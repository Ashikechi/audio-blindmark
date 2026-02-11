# cython: boundscheck = False
# cython: wraparound = False

from cpython.pycapsule cimport PyCapsule_GetPointer
from libc.stdint cimport *  # pyright: ignore[reportWildcardImportFromLibrary]

from numpy.random cimport bitgen

from .random import *


cpdef bytearray get_random_bytes(uint64_t size):
    cdef bitgen *rng = <bitgen*>PyCapsule_GetPointer(get_rng().capsule, <const char*>"BitGenerator")
    cdef bytearray r = bytearray(size)
    cdef uint8_t[::1] r_view = r # pyright: ignore[reportGeneralTypeIssues]
    for i in range(<uint64_t>len(r)): # pyright: ignore[reportGeneralTypeIssues]
        r_view[i] = <uint8_t>rng.next_uint32(rng.state)
    return r

cpdef bytearray encode(const uint8_t[::1] data):
    assert data is not None

    cdef bitgen *rng = <bitgen*>PyCapsule_GetPointer(get_rng().capsule, <const char*>"BitGenerator")
    cdef bytearray r = bytearray(data.shape[0] * 2)
    cdef uint8_t[::1] r_view = r # pyright: ignore[reportGeneralTypeIssues]
    cdef uint64_t i
    for i in range(<uint64_t>len(r)): # pyright: ignore[reportGeneralTypeIssues]
        if i & 1 == 0:
            r_view[i] = data[i // 2] >> 4 | <uint8_t>rng.next_uint32(rng.state) & ((<uint8_t>1 << 4) - 1) << 4 # pyright: ignore[reportGeneralTypeIssues]
        else:
            r_view[i] = data[i // 2] ^ (<uint8_t>rng.next_uint32(rng.state) & ((<uint8_t>1 << 4) - 1) << 4) # pyright: ignore[reportGeneralTypeIssues]
    return r

cpdef bytearray decode(const uint8_t[::1] data):
    assert data is not None
    assert data.shape[0] & 1 == 0

    cdef bytearray result = bytearray(data.shape[0] // 2)
    cdef int i
    for i in range(len(result)):
        result[i] = (data[i * 2] & (<uint8_t>1 << 4) - 1) << 4 | data[i * 2 | 1] & (<uint8_t>1 << 4) - 1 # pyright: ignore[reportGeneralTypeIssues]
    return result
