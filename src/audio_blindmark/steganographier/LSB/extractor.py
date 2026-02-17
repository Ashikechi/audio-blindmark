import numpy as np

from ...base import BaseExtractor
from .coder import decode


class LSBExtractor(BaseExtractor):
    def __init__(self, data_length: int) -> None:
        self.length = data_length
        super().__init__()

    def wave_length(self) -> int:
        return self.length << 3

    def data_length(self) -> int:
        return self.length

    def extract(self, wave: np.ndarray[tuple[int], np.dtype[np.float64]]) -> bytes:
        assert wave.shape[0] == self.wave_length()

        return decode(wave)
