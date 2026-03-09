from typing import Optional

import numpy as np
from scipy import fft

from ...base import BaseExtractor


class FFTExtractor(BaseExtractor):
    def __init__(self, wave_length: int = 4096, data_length: int = 64, center: Optional[int] = None, quantum: float = 0.2) -> None:
        self.__wave_length = wave_length
        self.__data_length = data_length
        self.quantum = quantum

        max_freq_idx = wave_length // 2
        bit_num = data_length << 3

        if center is None:
            center = max_freq_idx // 2
        if center - (bit_num - 1) // 2 < 1:
            center = (bit_num - 1) // 2 + 1
        elif center + bit_num // 2 + 1 >= max_freq_idx:
            center = max_freq_idx - (bit_num // 2 + 1) - 1

        self.points = list(range(center - (bit_num - 1) // 2, center + bit_num // 2 + 1))
        super().__init__()

    def wave_length(self) -> int:
        return self.__wave_length

    def data_length(self) -> int:
        return self.__data_length

    def extract(self, wave: np.ndarray) -> bytes:
        assert wave.shape[0] == self.wave_length()

        fft_result = fft.rfft(wave)
        phase = np.angle(fft_result)

        extracted_bits = []
        for pos in self.points:
            p = (phase[pos] + np.pi) % (2 * np.pi)
            q = int(p // self.quantum)
            extracted_bits.append(q % 2)

        return np.packbits(np.array(extracted_bits)).tobytes()
