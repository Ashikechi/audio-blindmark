from libc.stdint cimport * # pyright: ignore[reportWildcardImportFromLibrary]
from numpy.random cimport bitgen


cpdef bytearray get_random_bytes(uint64_t size)
cpdef bytearray encode(const uint8_t[::1] data)
cpdef bytearray decode(const uint8_t[::1] data)
