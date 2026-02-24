import sys
import os
from unittest.mock import MagicMock
import numpy as np
import wave
import pytest
from pytest_benchmark.fixture import BenchmarkFixture

def setup_final_mocks():
    m_blake3_obj = MagicMock()
    m_blake3_obj.digest.return_value = b'\x00' * 8
    m_blake3 = MagicMock(return_value=m_blake3_obj)
    sys.modules["blake3"] = MagicMock(blake3=m_blake3)

    m_creedsolo = MagicMock()
    m_creedsolo.RSCodec.return_value.encode.side_effect = lambda x: x
    m_creedsolo.RSCodec.return_value.decode.side_effect = lambda x: (x, None, None)
    sys.modules["creedsolo"] = m_creedsolo

    m_rt = MagicMock()
    m_rt.get_random_bytes.side_effect = lambda l: b'\x00' * l
    m_rt.encode.side_effect = lambda x: x + x
    m_rt.decode.side_effect = lambda x: x[:len(x)//2]
    sys.modules["audio_blindmark.utils.random_tools"] = m_rt

    m_xor = MagicMock()
    m_xor.xor.side_effect = lambda a, b: bytes([x ^ y for x, y in zip(a, b)])
    sys.modules["audio_blindmark.utils.bytearray_xor"] = m_xor

    m_random = MagicMock()
    m_random.get_rng.return_value = np.random.default_rng().bit_generator
    sys.modules["audio_blindmark.utils.random"] = m_random

setup_final_mocks()

from audio_blindmark import Encoder, Decoder, embed, extract
from audio_blindmark.steganographier.DCT import DCTEmbedder, DCTExtractor
from audio_blindmark.utils.random import seed

MEDIA_DIR = 'assets/audio/'
if not os.path.exists(MEDIA_DIR):
    os.makedirs(MEDIA_DIR)

KEY = b'Ema'
DATA = b'Sakuraba Ema daisuki' * 4
DATA_LENGTH = 64

def test_DCT():
    wav_path = os.path.join(MEDIA_DIR, 'test_short.wav')
    if os.path.exists(wav_path):
        os.remove(wav_path)

    with wave.open(wav_path, 'wb') as wf:
        wf.setparams((1, 2, 44100, 0, 'NONE', 'not compressed'))
        wf.writeframes(np.random.randint(-10000, 10000, 44100 * 5).astype(np.int16).tobytes())

    seed(42)
    embedder = DCTEmbedder(data_length=DATA_LENGTH, iterations=3)
    encoder = Encoder(KEY, ecc_length=0)

    output_path = os.path.join(MEDIA_DIR, 'output_dct.wav')

    with wave.open(wav_path, 'r') as raw_audio:
        with wave.open(output_path, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    extractor = DCTExtractor(data_length=DATA_LENGTH)
    decoder = Decoder(KEY, ecc_length=0)
    with wave.open(output_path, 'r') as audio:
        result = extract(audio, decoder, extractor)
        assert result == DATA

@pytest.mark.benchmark(group = 'DCT-embed')
def test_benchmark_DCT_embed(benchmark: BenchmarkFixture):
    def do():
        seed(42)
        embedder = DCTEmbedder(DATA_LENGTH, iterations=3)
        encoder = Encoder(KEY, ecc_length=0)
        with wave.open(MEDIA_DIR + 'test_short.wav', 'r') as raw_audio:
            with wave.open(MEDIA_DIR + 'output_dct.wav', 'w') as output_audio:
                output_audio.setparams(raw_audio.getparams())
                embed(raw_audio, DATA, output_audio, encoder, embedder)
    benchmark.pedantic(do, rounds=5)

@pytest.mark.benchmark(group = 'DCT-extract')
def test_benchmark_DCT_extract(benchmark: BenchmarkFixture):
    def do():
        extractor = DCTExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ecc_length=0)
        with wave.open(MEDIA_DIR + 'output_dct.wav', 'r') as audio:
            extract(audio, decoder, extractor)
    benchmark.pedantic(do, rounds=5)
