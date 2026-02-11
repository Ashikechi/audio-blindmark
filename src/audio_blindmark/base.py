from abc import ABC, abstractmethod

import numpy as np


class EmbedError(Exception):
    pass

class BaseEmbedder(ABC):
    @abstractmethod
    def wave_length(self) -> int:
        pass

    @abstractmethod
    def data_length(self) -> int:
        pass

    @abstractmethod
    def embed(self, wave: np.ndarray[tuple[int], np.dtype[np.float64]], data: bytes) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        pass

class ExtractError(Exception):
    pass

class BaseExtractor(ABC):
    @abstractmethod
    def wave_length(self) -> int:
        pass

    @abstractmethod
    def data_length(self) -> int:
        pass

    @abstractmethod
    def extract(self, wave: np.ndarray[tuple[int], np.dtype[np.float64]]) -> bytes:
        pass
