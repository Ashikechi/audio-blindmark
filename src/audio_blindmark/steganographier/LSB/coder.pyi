import numpy as np

def encode(wave: np.ndarray[tuple[int], np.dtype[np.float64]], data: bytes) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    ...

def decode(wave: np.ndarray[tuple[int], np.dtype[np.float64]]) -> bytearray:
    ...
