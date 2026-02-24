import numpy as np
from ...base import BaseEmbedder, EmbedError
from .coder import encode

class DCTEmbedder(BaseEmbedder):
    def __init__(self, data_length: int = 8, quantum: float = 1 / 4096, center: int = None, samples_per_byte: int = 1024, iterations: int = 3) -> None:
        self.length = data_length
        self.quantum = quantum
        self.center = center
        self.samples_per_byte = samples_per_byte
        self.iterations = iterations
        super().__init__()

    def wave_length(self) -> int:
        return self.length * self.samples_per_byte

    def data_length(self) -> int:
        return self.length

    def embed(self, wave: np.ndarray, data: bytes) -> np.ndarray:
        if wave.shape[0] < len(data) * 8:
            raise EmbedError(f"Segment length {wave.shape[0]} is too small to carry {len(data) * 8} bits")

        if np.sqrt(np.dot(wave, wave) / wave.shape[0]) <= 256:
            raise EmbedError("Segment energy too low")

        actual_center = self.center if self.center is not None else wave.shape[0] // 2

        try:
            return encode(wave, data, self.quantum, actual_center, self.iterations)
        except ValueError as e:
            raise EmbedError(str(e))
