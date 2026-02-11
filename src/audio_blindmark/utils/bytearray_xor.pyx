# cython: boundscheck = False
# cython: wraparound = False

from libc.stdint cimport *  # pyright: ignore[reportWildcardImportFromLibrary]


cpdef bytearray xor(const uint8_t[::1] a, const uint8_t[::1] b):
    assert a is not None and b is not None
    cdef uint64_t l = <uint64_t>(b.shape[0] if a.shape[0] > b.shape[0] else a.shape[0])
    cdef bytearray r = bytearray(l)
    cdef uint8_t[::1] r_view = r # pyright: ignore[reportGeneralTypeIssues]
    cdef uint64_t i
    for i in range(l): # pyright: ignore[reportGeneralTypeIssues]
        r_view[i] = <uint8_t>a[i] ^ <uint8_t>b[i] # pyright: ignore[reportGeneralTypeIssues]
    return r
