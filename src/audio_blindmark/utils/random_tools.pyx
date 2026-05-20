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
