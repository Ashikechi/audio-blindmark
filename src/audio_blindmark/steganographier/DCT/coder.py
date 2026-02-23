import numpy as np
import scipy.fft as fft

def encode(wave: np.ndarray, data: bytes) -> np.ndarray:
    quantum = 1 / 4096
    center = len(wave) // 2
    points = [i for i in range(center - 32, center + 32)]

    bits = []
    for byte in data:
        for i in range(8):
            bits.append((byte >> i) & 1)

    dct_result = np.array(fft.dct(wave, norm='ortho'))
    norm = np.linalg.norm(dct_result)

    dct_result /= norm
    now_norm = 1.0

    for pos, bit in zip(points, bits):
        current_val = dct_result[pos]
        if current_val // quantum % 2 == bit:
            new = (current_val // quantum + 0.5) * quantum
        else:
            new_neg = (current_val // quantum - 0.5) * quantum
            new_pos = (current_val // quantum + 1.5) * quantum
            if abs((now_norm - current_val**2 + new_neg**2) - 1) > \
               abs((now_norm - current_val**2 + new_pos**2) - 1):
                new = new_pos
            else:
                new = new_neg

        now_norm = now_norm - current_val**2 + new**2
        dct_result[pos] = new

    dct_result *= norm
    return np.array(fft.idct(dct_result, norm='ortho'))

def decode(wave: np.ndarray) -> bytes:
    quantum = 1 / 4096
    center = len(wave) // 2
    points = [i for i in range(center - 32, center + 32)]

    dct_result = np.array(fft.dct(wave, norm='ortho'))
    norm = np.linalg.norm(dct_result)
    dct_result /= norm

    extracted_bits = [int(dct_result[pos] // quantum % 2) for pos in points]

    res = bytearray()
    for i in range(0, len(extracted_bits), 8):
        byte = 0
        for j in range(8):
            byte |= (extracted_bits[i + j] << j)
        res.append(byte)
    return bytes(res)
