# cython: boundscheck = False
# cython: wraparound = False

from cpython.pycapsule cimport PyCapsule_GetPointer
from libc.stdint cimport *  # pylint: disable=W0401
from numpy.random cimport bitgen

from .random import *


cpdef bytearray get_random_bytes(uint64_t size):
    cdef bitgen *rng = <bitgen*>PyCapsule_GetPointer(get_rng().capsule, <const char*>"BitGenerator")
    cdef bytearray r = bytearray(size)
    cdef uint8_t[::1] r_view = r
    for i in range(<uint64_t>len(r)):
        r_view[i] = rng.next_uint32(rng.state)
    return r

cpdef bytearray encode(const uint8_t[::1] data):
    assert data is not None

    cdef bitgen *rng = <bitgen*>PyCapsule_GetPointer(get_rng().capsule, <const char*>"BitGenerator")
    cdef bytearray r = bytearray(data.shape[0] * 2)
    cdef uint8_t[::1] r_view = r
    cdef uint64_t i
    for i in range(<uint64_t>len(r)):
        if i & 1 == 0:
            r_view[i] = data[i // 2] >> 4 | <uint8_t>rng.next_uint32(rng.state) & ((1 << 4) - 1) << 4
        else:
            r_view[i] = data[i // 2] ^ (<uint8_t>rng.next_uint32(rng.state) & ((1 << 4) - 1) << 4)
    return r

cpdef bytearray decode(const uint8_t[::1] data):
    assert data is not None
    assert data.shape[0] & 1 == 0

    cdef bytearray result = bytearray(data.shape[0] // 2)
    cdef int i
    for i in range(len(result)):
        result[i] = (data[i * 2] & (1 << 4) - 1) << 4 | data[i * 2 | 1] & (1 << 4) - 1
    return result
