import numpy as np
from ...base import BaseExtractor, ExtractError
from .coder import decode

class DCTExtractor(BaseExtractor):
    def __init__(self, data_length: int = 8, quantum: float = 1/4096, center: int = None, samples_per_byte: int = 1024) -> None:
        self.length = data_length
        self.quantum = quantum
        self.center = center
        self.samples_per_byte = samples_per_byte
        super().__init__()

    def wave_length(self) -> int:
        return self.length * self.samples_per_byte

    def data_length(self) -> int:
        return self.length

    def extract(self, wave: np.ndarray) -> bytes:
        actual_center = self.center if self.center is not None else wave.shape[0] // 2

        if np.sqrt(np.dot(wave, wave) / wave.shape[0]) <= 256:
            raise ExtractError("No signal detected in this segment")

        return decode(wave, self.length, self.quantum, actual_center)
