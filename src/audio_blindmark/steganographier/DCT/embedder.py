import numpy as np
from ...base import BaseEmbedder, EmbedError
from .coder import encode

class DCTEmbedder(BaseEmbedder):
    def __init__(self, data_length: int = 8) -> None:
        """
        :param data_length: 每次处理嵌入的字节数。
        DCT 算法设定为每 1024 个采样点存 1 字节。
        默认 8 字节对应 8192 个采样点。
        """
        self.length = data_length
        super().__init__()

    def wave_length(self) -> int:
        return self.length * 1024

    def data_length(self) -> int:
        return self.length

    def embed(self, wave: np.ndarray, data: bytes) -> np.ndarray:
        assert wave.shape[0] == self.wave_length() and len(data) == self.data_length()

        if np.sqrt(np.dot(wave, wave) / wave.shape[0]) <= 256:
            raise EmbedError("Segment energy too low")

        return encode(wave, data)
