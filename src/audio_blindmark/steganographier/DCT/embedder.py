from typing import Optional

import numpy as np
from bitarray import bitarray
from scipy import fft

from ...base import BaseEmbedder, EmbedError


class DCTEmbedder(BaseEmbedder):
    def __init__(self, wave_length: int=4096, data_length: int=64, center: Optional[int]=None, quantum: float=1 / 1024, min_average_energy: float=256, iterations: int=3) -> None:
        self.__wave_length = wave_length
        self.__data_length = data_length
        self.quantum = quantum
        self.min_average_energy = min_average_energy
        self.iterations = iterations

        bit_num = data_length << 3
        if data_length << 3 > wave_length:
            raise ValueError(f'Can not embed {bit_num} bits in {wave_length} frames')

        if center is None:
            center = (wave_length - 1) // 2
        if center - (bit_num - 1) // 2 < 0:
            center = (bit_num - 1) // 2
        elif center + bit_num // 2 >= wave_length:
            center = wave_length - (bit_num // 2) - 1

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

        dct_result = fft.dct(wave, norm='ortho')
        average_energy = np.dot(dct_result, dct_result) / wave.shape[0]

        for _ in range(self.iterations):
            factor = np.sqrt(average_energy)
            for pos, bit in zip(self.points, bits):
                average_energy -= dct_result[pos] ** 2 / wave.shape[0]

                q = (dct_result[pos] / factor) // (self.quantum / 2)
                if q % 4 // 2 == bit:
                    q = q // 2 * 2 + 1
                else:
                    if q % 2 == 0:
                        q -= 1
                    else:
                        q += 2
                dct_result[pos] = q * (self.quantum / 2) * factor

                average_energy += dct_result[pos] ** 2 / wave.shape[0]

        return fft.idct(dct_result, norm='ortho')
