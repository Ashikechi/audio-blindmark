from typing import Optional

import numpy as np
from scipy import fft

from ...base import BaseExtractor


class DCTExtractor(BaseExtractor):
    def __init__(self, wave_length: int=4096, data_length: int=64, center: Optional[int]=None, quantum: float=1 / 4096, iterations: int=3) -> None:
        self.__wave_length = wave_length
        self.__data_length = data_length
        self.quantum = quantum
        self.iterations = iterations

        bit_num = data_length << 3
        if data_length << 3 > wave_length:
            raise ValueError(f'Can not extract {bit_num} bits from {wave_length} frames')

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

    def extract(self, wave: np.ndarray) -> bytes:
        assert wave.shape[0] == self.wave_length()

        dct_result = fft.dct(wave, norm='ortho')
        norm = np.linalg.norm(dct_result)
        dct_result /= norm

        return np.packbits((dct_result[self.points] // self.quantum).astype(np.int64) % 2).tobytes()
