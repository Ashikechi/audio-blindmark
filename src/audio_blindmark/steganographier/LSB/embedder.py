import numpy as np

from ...base import BaseEmbedder, EmbedError
from .coder import encode


class LSBEmbedder(BaseEmbedder):
    def __init__(self, data_length: int) -> None:
        self.length = data_length
        super().__init__()

    def wave_length(self) -> int:
        return self.length << 3

    def data_length(self) -> int:
        return self.length

    def embed(self, wave: np.ndarray[tuple[int], np.dtype[np.float64]], data: bytes) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        assert wave.shape[0] == self.wave_length() and len(data) == self.data_length()

        if np.sqrt(np.dot(wave, wave) / wave.shape[0]) <= 256:
            raise EmbedError
        return encode(wave, data)
