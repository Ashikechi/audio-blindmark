import numpy as np
import scipy.fft as fft
from bitarray import bitarray

def _encode_step(wave: np.ndarray, data_bits: bitarray, quantum: float, center: int) -> np.ndarray:
    num_bits = len(data_bits)
    half_width = num_bits // 2
    points = range(center - half_width, center + half_width)

    dct_result = np.array(fft.dct(wave, norm='ortho'))
    norm = np.linalg.norm(dct_result)

    dct_result /= norm
    now_norm = 1.0

    for pos, bit in zip(points, data_bits):
        current_val = dct_result[pos]

        if current_val // quantum % 2 == bit:
            new = (current_val // quantum + 0.5) * quantum
        else:
            new_neg = (current_val // quantum - 0.5) * quantum
            new_pos = (current_val // quantum + 1.5) * quantum

            if abs((now_norm - current_val ** 2 + new_neg ** 2) - 1) > \
               abs((now_norm - current_val ** 2 + new_pos ** 2) - 1):
                new = new_pos
            else:
                new = new_neg

        now_norm = now_norm - current_val ** 2 + new ** 2
        dct_result[pos] = new

    dct_result *= norm
    return np.array(fft.idct(dct_result, norm='ortho'))


def encode(wave: np.ndarray, data: bytes, quantum: float, center: int, iterations: int = 3) -> np.ndarray:
    num_bits = len(data) * 8
    half_width = num_bits // 2

    if center - half_width < 0 or center + half_width > len(wave):
        raise ValueError(f"Wave length {len(wave)} is too short for {len(data)} bytes at center {center}")

    data_bits = bitarray(endian='little')
    data_bits.frombytes(data)

    result = wave
    for _ in range(iterations):
        result = _encode_step(result, data_bits, quantum, center)
    return result

def decode(wave: np.ndarray, data_len: int, quantum: float, center: int) -> bytes:
    num_bits = data_len * 8
    half_width = num_bits // 2
    points = range(center - half_width, center + half_width)

    dct_result = np.array(fft.dct(wave, norm='ortho'))
    norm = np.linalg.norm(dct_result)
    dct_result /= norm

    extracted_bits = bitarray([int(dct_result[pos] // quantum % 2) for pos in points], endian='little')
    return extracted_bits.tobytes()
