import base64
from typing import Optional

from blake3 import blake3
from Crypto.Cipher import AES
from Crypto.Hash import SHA512
from Crypto.Protocol.KDF import PBKDF2
from Crypto.Random import get_random_bytes

from creedsolo import ReedSolomonError, RSCodec

from .utils import random_tools
from .utils.bytearray_xor import xor

DERIVATION_EPOCH = 1 << 20
DERIVATION_HASH = SHA512
MAX_DATA_LENGTH = 1 << 20


class ImpossibleOutputLengthError(Exception):
    pass

class Encoder:
    def __init__(self, password: bytes, ecc_length: int, digest_length = 4, salt_length = 4, nonce_length = 4) -> None:
        self.password = password
        self.ecc_length = ecc_length
        self.digest_length = digest_length
        self.salt_length = salt_length
        self.nonce_length = nonce_length

        self.rsc = RSCodec(ecc_length)

        key_over = PBKDF2(base64.b64encode(password).decode(), b'', 24, DERIVATION_EPOCH, hmac_hash_module = DERIVATION_HASH)
        cipher_over = AES.new(key_over, AES.MODE_CTR, nonce = b'')
        self.key_stream = cipher_over.encrypt(b'\x00' * MAX_DATA_LENGTH)

        self.key_cache: Optional[tuple[bytes, bytes]] = None

    def data_length(self, output_length: int) -> int:
        if 0 < output_length % 255 <= self.ecc_length:
            raise ImpossibleOutputLengthError
        return output_length - self.ecc_length * ((output_length - 1) // 255 + 1) - 4 - self.salt_length * 2 - self.nonce_length - self.digest_length

    def encode(self, seq: int, data: bytes) -> bytes:
        if self.key_cache:
            salt = self.key_cache[0]
        else:
            salt = get_random_bytes(self.salt_length)
            self.key_cache = (salt, PBKDF2(base64.b64encode(self.password).decode(), salt, 24, DERIVATION_EPOCH, hmac_hash_module = DERIVATION_HASH))

        nonce = get_random_bytes(self.nonce_length)

        cipher = AES.new(self.key_cache[1], AES.MODE_CTR, nonce = nonce)
        cipher_text = cipher.encrypt(seq.to_bytes(4, 'little') + data)
        cipher_text_with_meta = bytes(random_tools.encode(salt)) + nonce + cipher_text
        text_with_ecc: bytes = bytes(iter(self.rsc.encode(iter(blake3(cipher_text_with_meta).digest(self.digest_length) + cipher_text_with_meta))))

        assert len(text_with_ecc) < MAX_DATA_LENGTH
        return bytes(xor(text_with_ecc, self.key_stream))

class DecodeError(Exception):
    pass

class Decoder:
    def __init__(self, password: bytes, ecc_length: int, digest_length = 4, salt_length = 4, nonce_length = 4) -> None:
        self.password = password
        self.ecc_length = ecc_length
        self.digest_length = digest_length
        self.salt_length = salt_length
        self.nonce_length = nonce_length

        self.rsc = RSCodec(ecc_length)

        key_over = PBKDF2(base64.b64encode(password).decode(), b'', 24, DERIVATION_EPOCH, hmac_hash_module = DERIVATION_HASH)
        cipher_over = AES.new(key_over, AES.MODE_CTR, nonce = b'')
        self.key_stream = cipher_over.encrypt(b'\x00' * MAX_DATA_LENGTH)

        self.key_cache: dict[bytes, bytes] = {}

    def data_length(self, output_length: int) -> int:
        if 0 < output_length % 255 <= self.ecc_length:
            raise ImpossibleOutputLengthError
        return output_length - self.ecc_length * ((output_length - 1) // 255 + 1) - 4 - self.salt_length * 2 - self.nonce_length - self.digest_length

    def decode(self, data: bytes) -> tuple[int, bytes]:
        assert len(data) < MAX_DATA_LENGTH
        try:
            cipher_text_with_meta_with_digest: bytes = bytes(iter(self.rsc.decode(iter(xor(data, self.key_stream)))[0]))
        except ReedSolomonError as e:
            raise DecodeError(e) from e

        try:
            digest = cipher_text_with_meta_with_digest[:self.digest_length]
            cipher_text_with_meta = cipher_text_with_meta_with_digest[self.digest_length:]
            obfuscated_salt = cipher_text_with_meta[:self.salt_length * 2]
            nonce = cipher_text_with_meta[self.salt_length * 2:self.salt_length * 2 + self.nonce_length]
            cipher_text = cipher_text_with_meta[self.salt_length * 2 + self.nonce_length:]
        except IndexError as e:
            raise DecodeError(e) from e

        if blake3(cipher_text_with_meta).digest(self.digest_length) != digest:
            raise DecodeError('Unmatched digest')

        salt = bytes(random_tools.decode(obfuscated_salt))
        if salt in self.key_cache:
            key = self.key_cache[salt]
        else:
            self.key_cache[salt] = key = PBKDF2(base64.b64encode(self.password).decode(), salt, 24, DERIVATION_EPOCH, hmac_hash_module = DERIVATION_HASH)

        cipher = AES.new(key, AES.MODE_CTR, nonce = nonce)
        text = cipher.decrypt(cipher_text)
        try:
            return (int.from_bytes(text[:4], 'little'), text[4:])
        except IndexError as e:
            raise DecodeError(e) from e
