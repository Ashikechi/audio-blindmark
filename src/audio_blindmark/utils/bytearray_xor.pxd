from libc.stdint cimport * # pyright: ignore[reportWildcardImportFromLibrary]


cpdef bytearray xor(const uint8_t[::1] a, const uint8_t[::1] b)
