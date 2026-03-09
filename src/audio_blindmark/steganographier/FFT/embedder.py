from typing import Optional

import numpy as np
from bitarray import bitarray
from scipy import fft

from ...base import BaseEmbedder, EmbedError


class FFTEmbedder(BaseEmbedder):
    def __init__(self, wave_length: int=4096, data_length: int=64, center: Optional[int]=None, quantum: float=0.2, min_average_energy: float=256, iterations: int=3) -> None:
        self.__wave_length = wave_length
        self.__data_length = data_length
        self.quantum = quantum
        self.min_average_energy = min_average_energy
        self.iterations = iterations

        max_freq_idx = wave_length // 2
        bit_num = data_length << 3

        if bit_num > max_freq_idx - 2:
            raise ValueError(f'Can not embed {bit_num} bits in {wave_length} frames')

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

    def embed(self, wave: np.ndarray[tuple[int], np.dtype[np.float64]], data: bytes) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        assert wave.shape[0] == self.wave_length() and len(data) == self.data_length()

        if np.sqrt(np.dot(wave, wave) / wave.shape[0]) <= self.min_average_energy:
            raise EmbedError("Segment energy too low")

        bits = bitarray()
        bits.frombytes(data)

        current_wave = wave.copy()

        for _ in range(self.iterations):
            fft_result = fft.rfft(wave)
            magnitude = np.abs(fft_result)
            phase = np.angle(fft_result)

            for pos, bit in zip(self.points, bits):
                p = phase[pos] + np.pi
                q = p // self.quantum

                if q % 2 == bit:
                    new_p = (q + 0.5) * self.quantum
                else:
                    new_q = q + 1 if (p % self.quantum) > (self.quantum / 2) else q - 1
                    new_p = (new_q + 0.5) * self.quantum

                phase[pos] = (new_p % (2 * np.pi)) - np.pi

            new_fft_result = magnitude * np.exp(1j * phase)
            current_wave = fft.irfft(new_fft_result, n=self.__wave_length)

        return current_wave
