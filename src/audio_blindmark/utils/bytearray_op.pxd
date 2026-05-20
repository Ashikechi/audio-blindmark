from libc.stdint cimport *


cpdef bytearray xor(const uint8_t[::1] a, const uint8_t[::1] b)
cpdef bytearray add(const uint8_t[::1] a, const uint8_t[::1] b)
cpdef bytearray sub(const uint8_t[::1] a, const uint8_t[::1] b)
