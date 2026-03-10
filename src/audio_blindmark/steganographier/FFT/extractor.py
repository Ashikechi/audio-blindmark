from typing import Optional

import numpy as np
from scipy import fft

from ...base import BaseExtractor


class FFTExtractor(BaseExtractor):
    def __init__(self, wave_length: int = 4096, data_length: int = 64, center: Optional[int] = None, quantum: float = 0.8) -> None:
        self.__wave_length = wave_length
        self.__data_length = data_length
        self.quantum = quantum

        bit_num = data_length << 3
        max_bit_num = (wave_length - 1) // 2

        if bit_num > max_bit_num:
            raise ValueError(f'Can not extract {bit_num} bits from {wave_length} frames')

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

    def extract(self, wave: np.ndarray) -> bytes:
        assert wave.shape[0] == self.wave_length()

        fft_result = fft.rfft(wave, norm='ortho')
        angle = np.angle(fft_result)

        return np.packbits((angle[self.points] // self.quantum).astype(np.int64) % 2).tobytes()
