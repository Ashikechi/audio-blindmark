from libc.math cimport round
from libc.stdint cimport *  # pylint: disable=W0401

import numpy as np

from ...utils.random_tools import get_random_bytes


def encode(wave: np.ndarray[tuple[int], np.dtype[np.float64]], data: bytes) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    assert <uint64_t>wave.shape[0] == <uint64_t>len(data) << 3

    cdef const double[::1] wave_view = wave
    cdef const uint8_t[::1] data_view = data

    new_wave = np.empty(wave.shape[0])
    cdef double[::1] new_wave_view = new_wave

    cdef uint8_t[::1] random_bytes = get_random_bytes(data_view.shape[0])
    cdef uint64_t i
    cdef int64_t tmp
    for i in range(<uint64_t>new_wave_view.shape[0]):
        tmp = <int64_t>round(wave_view[i])
        if tmp & 1 == data_view[i >> 3] >> (i & ((1 << 3) - 1)) & 1:
            new_wave_view[i] = tmp
        else:
            new_wave_view[i] = tmp + <int64_t>(random_bytes[i >> 3] >> (i & ((1 << 3) - 1)) & 1) * 2 - 1
    return new_wave

def decode(wave: np.ndarray[tuple[int], np.dtype[np.float64]]) -> bytearray:
    assert wave.shape[0] & ((<uint64_t>1 << 3) - 1) == 0

    cdef const double[::1] wave_view = wave

    cdef bytearray data = bytearray(wave.shape[0] >> 3)
    cdef uint8_t[::1] data_view = data

    cdef uint64_t i
    for i in range(<uint64_t>wave_view.shape[0]):
        if <int64_t>round(wave_view[i]) & 1 == 1:
            data_view[i >> 3] |= 1 << (i & ((1 << 3) - 1))
    return data
