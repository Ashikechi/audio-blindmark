# cython: boundscheck = False
# cython: wraparound = False

from libc.stdint cimport *  # pylint: disable=W0401


cpdef bytearray xor(const uint8_t[::1] a, const uint8_t[::1] b):
    assert a is not None and b is not None
    cdef uint64_t l = b.shape[0] if a.shape[0] > b.shape[0] else a.shape[0]
    cdef bytearray r = bytearray(l)
    cdef uint8_t[::1] r_view = r
    cdef uint64_t i
    for i in range(l):
        r_view[i] = a[i] ^ b[i]
    return r
