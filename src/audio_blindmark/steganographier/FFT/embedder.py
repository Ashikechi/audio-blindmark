from typing import Optional

import numpy as np
from bitarray import bitarray
from scipy import fft

from ...base import BaseEmbedder, EmbedError


class FFTEmbedder(BaseEmbedder):
    def __init__(self, wave_length: int=4096, data_length: int=64, center: Optional[int]=None, quantum: float=0.8, min_average_energy: float=256) -> None:
        self.__wave_length = wave_length
        self.__data_length = data_length
        self.quantum = quantum
        self.min_average_energy = min_average_energy

        bit_num = data_length << 3
        max_bit_num = (wave_length - 1) // 2

        if bit_num > max_bit_num:
            raise ValueError(f'Can not embed {bit_num} bits in {wave_length} frames')

        if center is None:
            center = (max_bit_num + 1) // 2
        if center - (bit_num - 1) // 2 < 1:
            center = (bit_num - 1) // 2 + 1
        elif center + bit_num // 2 > max_bit_num:
            center = max_bit_num - (bit_num // 2)

        self.points = list(range(center - (bit_num - 1) // 2, center + bit_num // 2 + 1))
        super().__init__()

    def wave_length(self) -> int:
        return self.__wave_length

    def data_length(self) -> int:
        return self.__data_length

    def embed(self, wave: np.ndarray[tuple[int], np.dtype[np.float64]], data: bytes) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        assert wave.shape[0] == self.wave_length() and len(data) == self.data_length()

        if np.sqrt(np.dot(wave, wave) / wave.shape[0]) <= self.min_average_energy:
            raise EmbedError("Segment energy too low")

        bits = bitarray(data)

        fft_result = fft.rfft(wave, norm='ortho')

        for pos, bit in zip(self.points, bits):
            magnitude = np.abs(fft_result[pos])
            angle_quantized = np.angle(fft_result[pos]) // (self.quantum / 2)

            if angle_quantized % 4 // 2 == bit:
                angle_quantized = angle_quantized // 2 * 2 + 1
            else:
                if angle_quantized % 2 == 0:
                    angle_quantized -= 1
                else:
                    angle_quantized += 2

            fft_result[pos] = magnitude * np.exp(1j * angle_quantized * (self.quantum / 2))

        return fft.irfft(fft_result, wave.shape[0], norm='ortho')
