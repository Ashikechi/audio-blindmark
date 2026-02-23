import numpy as np
from ...base import BaseExtractor, ExtractError
from .coder import decode

class DCTExtractor(BaseExtractor):
    def __init__(self, data_length: int = 8) -> None:
        self.length = data_length
        super().__init__()

    def wave_length(self) -> int:
        return self.length * 1024

    def data_length(self) -> int:
        return self.length

    def extract(self, wave: np.ndarray) -> bytes:
        assert wave.shape[0] == self.wave_length()

        if np.sqrt(np.dot(wave, wave) / wave.shape[0]) <= 256:
            raise ExtractError("No signal detected in this segment")

        return decode(wave)
