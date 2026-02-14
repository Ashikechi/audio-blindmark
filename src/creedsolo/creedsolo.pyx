# cython: language_level = 3str

# -*- coding: utf-8 -*-

# Copyright (c) 2012-2015 Tomer Filiba <tomerfiliba@gmail.com>
# Copyright (c) 2015 rotorgit
# Copyright (c) 2015-2023 Stephen Larroque <LRQ3000@gmail.com>
# Copyright (c) 2012-2015 Night Space <NightSpaceC@outlook.com>

'''
Reed Solomon
============

A pure-python `universal errors-and-erasures Reed-Solomon Codec <http://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction>`_
, based on the wonderful tutorial at
`wikiversity <http://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders>`_,
written by "Bobmath" and "LRQ3000".

The code of wikiversity is here consolidated into a nice API with exceptions handling.
The algorithm can correct up to 2*e+v <= nsym, where e is the number of errors,
v the number of erasures and nsym = n-k = the number of ECC (error correction code) symbols.
This means that you can either correct exactly floor(nsym/2) errors, or nsym erasures
(errors where you know the position), and a combination of both errors and erasures.
The code should work on pretty much any reasonable version of python (2.4-3.5),
but I'm only testing on 2.7 - 3.4.

.. note::
   The codec is universal, meaning that it can decode any message encoded by another RS encoder
   as long as you provide the correct parameters.
   Note however that if you use higher fields (ie, bigger c_exp), the algorithms will be slower, first because
   we cannot then use the optimized bytearray() structure but only pa.array('i', ...), and also because
   Reed-Solomon's complexity is quadratic (both in encoding and decoding), so this means that the longer
   your messages, the longer it will take to encode/decode (quadratically!).

   The algorithm itself can handle messages up to (2^c_exp)-1 symbols, including the ECC symbols,
   and each symbol can have a value of up to (2^c_exp)-1 (indeed, both the message length and the maximum
   value for one character is constrained by the same mathematical reason). By default, we use the field GF(2^8),
   which means that you are limited to values between 0 and 255 (perfect to represent a single hexadecimal
   symbol on computers, so you can encode any binary stream) and limited to messages+ecc of maximum
   length 255. However, you can "chunk" longer messages to fit them into the message length limit.
   The ``RSCodec`` class will automatically apply chunking, by splitting longer messages into chunks and
   encode/decode them separately; it shouldn't make a difference from an API perspective (ie, from your POV).

::

    # Initialization
    >>> from creedsolo import RSCodec
    >>> rsc = RSCodec(10)  # 10 ecc symbols

    # Encoding
    >>> rsc.encode([1, 2, 3, 4])
    array('Q', [1, 2, 3, 4, 44, 157, 28, 43, 61, 248, 104, 250, 152, 77])
    >>> bytes(iter(rsc.encode(iter(b'hello world'))))
    b'hello world\xed%T\xc4\xfd\xfd\x89\xf3\xa8\xaa'
    # Note that chunking is supported transparently to encode any string length.

    # Decoding (repairing)
    >>> bytes(iter(rsc.decode(iter(b'hello world\xed%T\xc4\xfd\xfd\x89\xf3\xa8\xaa'))[0]))
    b'hello world'
    >>> bytes(iter(rsc.decode(iter(b'heXlo worXd\xed%T\xc4\xfdX\x89\xf3\xa8\xaa'))[0]))     # 3 errors
    b'hello world'
    >>> bytes(iter(rsc.decode(iter(b'hXXXo worXd\xed%T\xc4\xfdX\x89\xf3\xa8\xaa'))[0]))     # 5 errors
    b'hello world'
    >>> bytes(iter(rsc.decode(iter(b'hXXXo worXd\xed%T\xc4\xfdXX\xf3\xa8\xaa'))[0]))        # 6 errors - fail
    Traceback (most recent call last):
      ...
    ReedSolomonError: Could not locate error

    >>> rsc = RSCodec(12)  # using 2 more ecc symbols (to correct max 6 errors or 12 erasures)
    >>> bytes(iter(rsc.encode(iter(b'hello world'))))
    b'hello world?Ay\xb2\xbc\xdc\x01q\xb9\xe3\xe2='
    >>> bytes(iter(rsc.decode(iter(b'hello worXXXXy\xb2XX\x01q\xb9\xe3\xe2='))[0]))         # 6 errors - ok
    b'hello world'
    >>> bytes(iter(rsc.decode(iter(b'helXXXXXXXXXXy\xb2XX\x01q\xb9\xe3\xe2='), erase_pos=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16])[0]))  # 12 erasures - OK
    b'hello world'

    # Checking
    >> rsc.check(iter(b'hello worXXXXy\xb2XX\x01q\xb9\xe3\xe2='))
    array('B', [0])
    >> rmes, rmesecc, _ = rsc.decode(iter(b'hello worXXXXy\xb2XX\x01q\xb9\xe3\xe2='))
    >> rsc.check(rmesecc)
    array('B', [1])

    # To use longer chunks or bigger values than 255 (may be very slow)
    >> rsc = RSCodec(12, nsize = 4095)  # always use a power of 2 minus 1
    >> rsc = RSCodec(12, c_exp = 12)  # alternative way to set nsize=4095
    >> mes = b'a' * (4095-12)
    >> mesecc = rsc.encode(iter(mes))
    >> mesecc[2] = 1
    >> mesecc[-1] = 1
    >> rmes, rmesecc, _ = rsc.decode(mesecc)
    >> rsc.check(mesecc)
    array('B', [0])
    >> rsc.check(rmesecc)
    array('B', [1])

    If you want full control, you can skip the API and directly use the library as-is. Here's how:

    First you need to init the precomputed tables:
    >> import creedsolo.reedsolo as rs
    >> gf = rs.GaloisField(0x11d)
    Pro tip: if you get the error: ValueError: byte must be in range(256), please check that your prime polynomial is correct for your field.
    Pro tip2: by default, you can only encode messages of max length and max symbol value = 256. If you want to encode bigger messages,
    please use the following (where c_exp is the exponent of your Galois Field, eg, 12 = max length 2^12 = 4096):
    >> prim = rs.find_prime_polys(c_exp=12, fast_primes=True, single=True)[0]
    >> gf = rs.GaloisField(c_exp = 12, prim = prim)

    Let's define our RS message and ecc size:
    >> n = 255  # length of total message+ecc
    >> nsym = 12  # length of ecc
    >> mes = b'a' * (n - nsym)  # generate a sample message

    To optimize, you can precompute the generator polynomial:
    >> gen = rs.rs_generator_poly_all(gf, n)

    Then to encode:
    >> mesecc = rs.rs_encode_msg(gf, iter(mes), nsym, gen = gen[nsym])

    Let's tamper our message:
    >> mesecc[1] = 0

    To decode:
    >> rmes, recc, errata_pos = rs.rs_correct_msg(gf, mesecc, nsym, erase_pos = erase_pos)
    Note that both the message and the ecc are corrected (if possible of course).
    Pro tip: if you know a few erasures positions, you can specify them in a list `erase_pos` to double the repair power. But you can also just specify an empty list.

    If the decoding fails, it will normally automatically check and raise a ReedSolomonError exception that you can handle.
    However if you want to manually check if the repaired message is correct, you can do so:
    >> rs.rs_check(gf, rmes + recc, nsym)

    Read the sourcecode's comments for more info about how it works, and for the various parameters you can setup if
    you need to interface with other RS codecs.

'''

# TODO IMPORTANT: try to keep the same convention for the ordering of polynomials inside lists throughout the code and functions (because for now there are a lot of list reversing in order to make it work, you never know the order of a polynomial, ie, if the first coefficient is the major degree or the constant term...).

import array as pa
from typing import Annotated, Generator, Iterable, Optional, Sequence

cimport cpython.array as ca
from libc.stdint cimport *

################### INIT and stuff ###################

class ReedSolomonError(Exception):
    pass

cdef ca.array uint64_array_template = <ca.array>pa.array('Q')
cdef ca.array uint8_array_template = <ca.array>pa.array('B')


################### GALOIS FIELD ELEMENTS MATHS ###################
# General note: Galois Field maths essentially are all the standard math operations everybody learn in primary school,
# but with only integer AND they can wraparound (ie, we use a modulo), so that in practice this means that
# Galois Field math operations are bounded in a very specific range of values. This changes a lot how the maths are done,
# but not that much, so you can still wrap your head around if you are willing to spend some time.
cpdef ca.array rwh_primes1(uint64_t n):
    ''' Returns a list of primes < n '''
    cdef ca.array is_prime = ca.clone(uint8_array_template, n, False)
    cdef uint64_t i
    for i in range(2, <uint64_t>len(is_prime)):
        is_prime.data.as_uchars[i] = 1
    is_prime.data.as_uchars[0] = is_prime.data.as_uchars[1] = 0

    cdef ca.array primes = ca.clone(uint64_array_template, n, False)
    cdef uint64_t primes_size = 0
    cdef uint64_t j
    for i in range(2, n):
        if is_prime.data.as_uchars[i]:
            primes.data.as_ulonglongs[primes_size] = i
            primes_size += 1
        for j in range(<uint64_t>len(primes)):
            if i * primes.data.as_ulonglongs[j] >= n:
                break
            is_prime.data.as_uchars[i * primes.data.as_ulonglongs[j]] = 0
            if i % primes.data.as_ulonglongs[j] == 0:
                break
    ca.resize(primes, primes_size)
    return primes

# ### Define bitwise carry-less operations as inner functions ###
# cpdef uint64_t cl_mult(uint64_t x, uint64_t y):
#     '''Bitwise carry-less multiplication on integers'''
#     cdef uint64_t z = 0
#     cdef uint64_t i = 0
#     while (y >> i) > 0:
#         if y & (1 << i):
#             z ^= x << i
#         i += 1
#     return z

# cpdef uint64_t bit_length(uint64_t n):
#     '''Compute the position of the most significant bit (1) of an integer. Equivalent to int.bit_length()'''
#     cdef uint64_t bits = 0
#     while n >> bits:
#         bits += 1
#     return bits

# cpdef uint64_t cl_div(uint64_t dividend, uint64_t divisor):
#     '''Bitwise carry-less long division on integers and returns the remainder'''
#     # Compute the position of the most significant bit for each integers
#     cdef uint64_t dl1 = bit_length(dividend)
#     cdef uint64_t dl2 = bit_length(divisor)
#     # If the dividend is smaller than the divisor, just exit
#     if dl1 < dl2:
#         return dividend
#     # Else, align the most significant 1 of the divisor to the most significant 1 of the dividend (by shifting the divisor)
#     cdef uint64_t i
#     for i in range(dl1 - dl2, -1, -1):
#         # Check that the dividend is divisible (useless for the first iteration but important for the next ones)
#         if dividend & 1 << i + dl2 - 1:
#             # If divisible, then shift the divisor to align the most significant bits and XOR (carry-less substraction)
#             dividend ^= divisor << i
#     return dividend

# ### Main GF multiplication routine ###
# cpdef uint64_t gf_mult_noLUT_slow(uint64_t x, uint64_t y, uint64_t prim = 0):  # pylint: disable=C0103
#     '''Multiplication in Galois Fields on-the-fly without using a precomputed look-up table (and thus it's slower) by using the standard carry-less multiplication + modular reduction using an irreducible prime polynomial.'''
#     # Multiply the gf numbers
#     result = cl_mult(x, y)
#     # Then do a modular reduction (ie, remainder from the division) with an irreducible primitive polynomial so that it stays inside GF bounds
#     if prim != 0:
#         result = cl_div(result, prim)

#     return result

cpdef uint64_t gf_mult_noLUT(uint64_t x, uint64_t y, uint64_t prim = 0, uint64_t field_charac_full = 256, bool carryless = True):  # pylint: disable=C0103
    '''Galois Field integer multiplication on-the-fly without using a look-up table, using Russian Peasant Multiplication algorithm (faster than the standard multiplication + modular reduction). This is still slower than using a look-up table, but is the fastest alternative, and is often used in embedded circuits where storage space is limited (ie, no space for a look-up table).
    If prim is 0 and carryless=False, then the function produces the result for a standard integers multiplication (no carry-less arithmetics nor modular reduction).'''
    cdef uint64_t r = 0
    while y: # while y is above 0
        if y & 1 != 0:
            r = r ^ x if carryless else r + x  # y is odd, then add the corresponding x to r (the sum of all x's corresponding to odd y's will give the final product). Note that since we're in GF(2), the addition is in fact an XOR (very important because in GF(2) the multiplication and additions are carry-less, thus it changes the result!).
        y = y >> 1  # equivalent to y // 2
        x = x << 1  # equivalent to x * 2
        if prim > 0 and x & field_charac_full != 0:
            x = x ^ prim  # GF modulo: if x >= 256 then apply modular reduction using the primitive polynomial (we just substract, but since the primitive number can be above 256 then we directly XOR).

    return r

cpdef ca.array find_prime_polys(uint64_t generator = 2, uint64_t c_exp = 8, bool fast_primes = False, bool single = False):
    '''Compute the list of prime polynomials for the given generator and galois field characteristic exponent.'''
    # fast_primes will output less results but will be significantly faster.
    # single will output the first prime polynomial found, so if all you want is to just find one prime polynomial to generate the LUT for Reed-Solomon to work, then just use that.

    # A prime polynomial (necessarily irreducible) is necessary to reduce the multiplications in the Galois Field, so as to avoid overflows.
    # Why do we need a "prime polynomial"? Can't we just reduce modulo 255 (for GF(2^8) for example)? Because we need the values to be unique.
    # For example: if the generator (alpha) = 2 and c_exp = 8 (GF(2^8) == GF(256)), then the generated Galois Field (0, 1, α, α^1, α^2, ..., α^(p-1)) will be galois field it becomes 0, 1, 2, 4, 8, 16, etc. However, upon reaching 128, the next value will be doubled (ie, next power of 2), which will give 256. Then we must reduce, because we have overflowed above the maximum value of 255. But, if we modulo 255, this will generate 256 == 1. Then 2, 4, 8, 16, etc. giving us a repeating pattern of numbers. This is very bad, as it's then not anymore a bijection (ie, a non-zero value doesn't have a unique index). That's why we can't just modulo 255, but we need another number above 255, which is called the prime polynomial.
    # Why so much hassle? Because we are using precomputed look-up tables for multiplication: instead of multiplying a*b, we precompute alpha^a, alpha^b and alpha^(a+b), so that we can just use our lookup table at alpha^(a+b) and get our result. But just like in our original field we had 0,1,2,...,p-1 distinct unique values, in our "LUT" field using alpha we must have unique distinct values (we don't care that they are different from the original field as long as they are unique and distinct). That's why we need to avoid duplicated values, and to avoid duplicated values we need to use a prime irreducible polynomial.

    # Here is implemented a bruteforce approach to find all these prime polynomials, by generating every possible prime polynomials (ie, every integers between field_charac+1 and field_charac*2), and then we build the whole Galois Field, and we reject the candidate prime polynomial if it duplicates even one value or if it generates a value above field_charac (ie, cause an overflow).
    # Note that this algorithm is slow if the field is too big (above 12), because it's an exhaustive search algorithm. There are probabilistic approaches, and almost surely prime approaches, but there is no determistic polynomial time algorithm to find irreducible monic polynomials. More info can be found at: http://people.mpi-inf.mpg.de/~csaha/lectures/lec9.pdf
    # Another faster algorithm may be found at Adleman, Leonard M., and Hendrik W. Lenstra. "Finding irreducible polynomials over finite fields." Proceedings of the eighteenth annual ACM symposium on Theory of computing. ACM, 1986.

    # Prepare the finite field characteristic (2^p - 1), this also represent the maximum possible value in this field
    # we're in GF(2)
    cdef uint64_t field_charac = (<uint64_t>1 << c_exp) - 1
    cdef uint64_t field_charac_next = (<uint64_t>1 << c_exp + 1) - 1
    cdef ca.array prim_candidates
    cdef uint64_t prim_candidates_size = 0
    cdef ca.array primes
    cdef uint64_t i
    if fast_primes:
        primes = rwh_primes1(field_charac_next)  # generate maybe prime polynomials and check later if they really are irreducible
        prim_candidates = ca.clone(uint64_array_template, len(primes), False)
        for i in range(<uint64_t>len(primes)):  # filter out too small primes
            if primes.data.as_ulonglongs[i] > field_charac:
                prim_candidates.data.as_ulonglongs[prim_candidates_size] = primes.data.as_ulonglongs[i]
                prim_candidates_size += 1
    else:
        prim_candidates = ca.clone(uint64_array_template, <Py_ssize_t>(1 << c_exp - 1) - 1, False)
        for i in range(field_charac + 2, field_charac_next, 2):  # try each possible prime polynomial, but skip even numbers (because divisible by 2 so necessarily not irreducible)
            prim_candidates.data.as_ulonglongs[prim_candidates_size] = i
            prim_candidates_size += 1
    ca.resize(prim_candidates, prim_candidates_size)

    # Start of the main loop
    cdef ca.array correct_primes = ca.clone(uint64_array_template, len(prim_candidates), False)
    cdef uint64_t correct_primes_size = 0
    cdef ca.array seen
    cdef bool conflict
    cdef uint64_t _, x
    for i in range(<uint64_t>len(prim_candidates)):  # try potential candidates primitive irreducible polys
        seen = ca.clone(uint64_array_template, <Py_ssize_t>(field_charac + 1), True)  # memory variable to indicate if a value was already generated in the field (value at index x is set to 1) or not (set to 0 by default)
        conflict = False  # flag to know if there was at least one conflict

        # Second loop, build the whole Galois Field
        x = 1
        for _ in range(field_charac):
            # Compute the next value in the field (ie, the next power of alpha/generator)
            x = gf_mult_noLUT(x, generator, prim_candidates.data.as_ulonglongs[i], field_charac + 1)

            # Rejection criterion: if the value overflowed (above field_charac) or is a duplicate of a previously generated power of alpha, then we reject this polynomial (not prime)
            if x > field_charac or seen.data.as_ulonglongs[x] == 1:
                conflict = True
                break
            # Else we flag this value as seen (to maybe detect future duplicates), and we continue onto the next power of alpha
            seen.data.as_ulonglongs[x] = 1

        # End of the second loop: if there's no conflict (no overflow nor duplicated value), this is a prime polynomial!
        if not conflict:
            if single:
                ca.resize(correct_primes, 1)
                correct_primes.data.as_ulonglongs[0] = prim_candidates.data.as_ulonglongs[i]
                return correct_primes  # for API consistency, we always return an array, but here with a single value
            correct_primes.data.as_ulonglongs[correct_primes_size] = prim_candidates.data.as_ulonglongs[i]
            correct_primes_size += 1

    ca.resize(correct_primes, correct_primes_size)
    # Return the list of all prime polynomials
    return correct_primes  # you can use the following to print the hexadecimal representation of each prime polynomial: print [hex(i) for i in correct_primes]

cdef class GaloisField:
    def __init__(self, prim: int = 0x11d, generator: int = 2, c_exp: int = 8) -> None:
        '''Precompute the logarithm and anti-log tables for faster computation later, using the provided primitive polynomial.
        These tables are used for multiplication/division since addition/substraction are simple XOR operations inside GF of characteristic 2.
        The basic idea is quite simple: since b ** (log_b(x), log_b(y)) == x * y given any number b (the base or generator of the logarithm), then we can use any number b to precompute logarithm and anti-log (exponentiation) tables to use for multiplying two numbers x and y.
        That's why when we use a different base/generator number, the log and anti-log tables are drastically different, but the resulting computations are the same given any such tables.
        For more infos, see https://en.wikipedia.org/wiki/Finite_field_arithmetic#Implementation_tricks
        '''
        # generator is the generator number (the "increment" that will be used to walk through the field by multiplication, this must be a prime number). This is basically the base of the logarithm/anti-log tables. Also often noted "alpha" in academic books.
        # prim is the primitive/prime (binary) polynomial and must be irreducible (ie, it can't represented as the product of two smaller polynomials). It's a polynomial in the binary sense: each bit is a coefficient, but in fact it's an integer between field_charac+1 and field_charac*2, and not a list of gf values. The prime polynomial will be used to reduce the overflows back into the range of the Galois Field without duplicating values (all values should be unique). See the function find_prime_polys() and: http://research.swtch.com/field and http://www.pclviewer.com/rs2/galois.html
        # note that the choice of generator or prime polynomial doesn't matter very much: any two finite fields of size p^n have identical structure, even if they give the individual elements different names (ie, the coefficients of the codeword will be different, but the final result will be the same: you can always correct as many errors/erasures with any choice for those parameters). That's why it makes sense to refer to all the finite fields, and all decoders based on Reed-Solomon, of size p^n as one concept: GF(p^n). It can however impact sensibly the speed (because some parameters will generate sparser tables).
        # c_exp is the exponent for the field's characteristic GF(2^c_exp)

        # Init global tables
        cdef uint64_t field_charac = (1 << c_exp) - 1
        cdef ca.array gf_exp = ca.clone(uint64_array_template, <Py_ssize_t>(field_charac * 2), False) # anti-log (exponential) table. The first two elements will always be [GF256int(1), generator]
        cdef ca.array gf_log = ca.clone(uint64_array_template, <Py_ssize_t>(field_charac + 1), False) # log table, log[0] is impossible and thus unused

        # For each possible value in the galois field 2^8, we will pre-compute the logarithm and anti-logarithm (exponential) of this value
        # To do that, we generate the Galois Field F(2^p) by building a list starting with the element 0 followed by the (p-1) successive powers of the generator α : 1, α, α^1, α^2, ..., α^(p-1).
        cdef uint64_t x = 1
        cdef uint64_t i
        for i in range(field_charac):  # we could skip index 255 which is equal to index 0 because of modulo: g^255==g^0 but either way, this does not change the later outputs (ie, the ecc symbols will be the same either way)
            gf_exp.data.as_ulonglongs[i] = x  # compute anti-log for this value and store it in a table
            gf_log.data.as_ulonglongs[x] = i  # compute log at the same time
            x = gf_mult_noLUT(x, generator, prim, field_charac + 1)

            # If you use only generator==2 or a power of 2, you can use the following which is faster than gf_mult_noLUT():
            #x <<= 1  # multiply by 2 (change 1 by another number y to multiply by a power of 2^y)
            #if x & 0x100:  # similar to x >= 256, but a lot faster (because 0x100 == 256)
                #x ^= prim  # substract the primary polynomial to the current value (instead of 255, so that we get a unique set made of coprime numbers), this is the core of the tables generation

        # Optimization: double the size of the anti-log table so that we don't need to mod 255 to stay inside the bounds (because we will mainly use this table for the multiplication of two GF numbers, no more).
        for i in range(field_charac, field_charac * 2):
            gf_exp.data.as_ulonglongs[i] = gf_exp.data.as_ulonglongs[i - field_charac]
        self.field_charac, self.gf_exp, self.gf_log = field_charac, gf_exp, gf_log

    cpdef uint64_t gf_inverse(self, uint64_t x):
        '''Inverse of a galois field integer'''
        return self.gf_exp.data.as_ulonglongs[self.field_charac - self.gf_log.data.as_ulonglongs[x]]  # gf_inverse(x) == gf_div(1, x)

    cpdef uint64_t gf_mul(self, uint64_t x, uint64_t y):
        '''Multiply two galois field integers'''
        if x == 0 or y == 0:
            return 0
        return self.gf_exp.data.as_ulonglongs[(self.gf_log.data.as_ulonglongs[x] + self.gf_log.data.as_ulonglongs[y]) % self.field_charac]

    cpdef uint64_t gf_div(self, uint64_t x, uint64_t y):
        '''Divide x by y galois field integers'''
        if y == 0:
            raise ZeroDivisionError()
        if x == 0:
            return 0
        return self.gf_exp.data.as_ulonglongs[(self.gf_log.data.as_ulonglongs[x] + self.field_charac - self.gf_log.data.as_ulonglongs[y]) % self.field_charac]

    cpdef uint64_t gf_pow(self, uint64_t x, int64_t power):
        '''Power of x galois field integer'''
        return self.gf_exp.data.as_ulonglongs[(<int64_t>self.gf_log.data.as_ulonglongs[x] * power) % <int64_t>self.field_charac]

cpdef uint64_t gf_add(uint64_t x, uint64_t y):
    '''Add two galois field integers'''
    return x ^ y

cpdef uint64_t gf_sub(uint64_t x, uint64_t y):
    '''Subtract two galois field integers'''
    return x ^ y  # in binary galois field, substraction is just the same as addition (since we mod 2)

cpdef uint64_t gf_neg(uint64_t x):
    '''Negate one galois field integer (does nothing)'''
    return x


################### GALOIS FIELD POLYNOMIALS MATHS ###################

cpdef ca.array gf_poly_scale(GaloisField gf, const uint64_t[::1] p, uint64_t x):
    '''Scale a galois field polynomial with a factor x (an integer)'''
    assert p is not None
    #return pa.array('Q', (gf_mul(p[i], x) for i in range(len(p))))  # unoptimized one-liner
    cdef ca.array out = ca.clone(uint64_array_template, p.shape[0], False)
    cdef uint64_t i
    for i in range(<uint64_t>p.shape[0]):
        out.data.as_ulonglongs[i] = gf.gf_mul(p[i], x)
    return out

cpdef ca.array gf_poly_add(const uint64_t[::1] p, const uint64_t[::1] q):
    '''Add two galois field polynomials'''
    assert p is not None and q is not None
    cdef uint64_t p_len = p.shape[0]
    cdef uint64_t q_len = q.shape[0]
    cdef ca.array r
    cdef uint64_t[::1] r_view
    cdef uint64_t i
    if q_len > p_len:
        r = ca.clone(uint64_array_template, q_len, False)
        r_view = r
        r_view[:] = q
        for i in range(p_len):
            r.data.as_ulonglongs[r_view.shape[0] - p_len + i] ^= p[i]
    else:
        r = ca.clone(uint64_array_template, p_len, False)
        r_view = r
        r_view[:] = p
        for i in range(q_len):
            r.data.as_ulonglongs[r_view.shape[0] - q_len + i] ^= q[i]
    return r

cpdef ca.array gf_poly_mul(GaloisField gf, const uint64_t[::1] p, const uint64_t[::1] q):
    '''Multiply two polynomials, inside Galois Field (but the procedure is generic). Optimized function by precomputation of log.'''
    assert p is not None and q is not None
    # Pre-allocate the result array
    cdef ca.array r = ca.clone(uint64_array_template, <Py_ssize_t>(p.shape[0] + q.shape[0] - 1), True)
    # Precompute the logarithm of p
    cdef ca.array lp = ca.clone(uint64_array_template, p.shape[0], False)
    cdef uint64_t i
    for i in range(<uint64_t>p.shape[0]):
        lp.data.as_ulonglongs[i] = gf.gf_log.data.as_ulonglongs[p[i]]
    # Compute the polynomial multiplication (just like the outer product of two vectors, we multiply each coefficients of p with all coefficients of q)
    cdef uint64_t j
    for j in range(<uint64_t>q.shape[0]):
        qj = q[j]  # optimization: load the coefficient once
        if qj != 0:  # log(0) is undefined, we need to check that
            lq = gf.gf_log.data.as_ulonglongs[qj]  # Optimization: precache the logarithm of the current coefficient of q
            for i in range(<uint64_t>p.shape[0]):
                if p[i] != 0:  # log(0) is undefined, need to check that...
                    r.data.as_ulonglongs[i + j] ^= gf.gf_exp[lp.data.as_ulonglongs[i] + lq]  # equivalent to: r[i + j] = gf_add(r[i + j], gf_mul(p[i], q[j]))
    return r

# cpdef ca.array gf_poly_mul_simple(GaloisField gf, const uint64_t[::1] p, const uint64_t[::1] q): # simple equivalent way of multiplying two polynomials without precomputation, but thus it's slower
#     '''Multiply two polynomials, inside Galois Field'''
#     assert p is not None and q is not None
#     # Pre-allocate the result array
#     cdef ca.array r = ca.clone(uint64_array_template, <Py_ssize_t>(p.shape[0] + q.shape[0] - 1), False)
#     # Compute the polynomial multiplication (just like the outer product of two vectors, we multiply each coefficients of p with all coefficients of q)
#     cdef uint64_t i, j
#     for i in range(<uint64_t>p.shape[0]):
#         for j in range(<uint64_t>q.shape[0]):
#             r.data.as_ulonglongs[i + j] ^= gf.gf_mul(p[i], q[j])  # equivalent to: r[i + j] = gf_add(r[i + j], gf_mul(p[i], q[j])) -- you can see it's your usual polynomial multiplication
#     return r

cpdef ca.array gf_poly_neg(const uint64_t[::1] poly):
    '''Returns the polynomial with all coefficients negated. In GF(2^p), negation does not change the coefficient, so we return the polynomial as-is.'''
    assert poly is not None
    cdef ca.array r = ca.clone(uint64_array_template, poly.shape[0], False)
    cdef uint64_t[::1] r_view = r
    r_view[...] = poly
    return r

cpdef tuple gf_poly_div(GaloisField gf, const uint64_t[::1] dividend, const uint64_t[::1] divisor):
    '''Fast polynomial division by using Extended Synthetic Division and optimized for GF(2^p) computations (doesn't work with standard polynomials outside of this galois field).'''
    assert dividend is not None and divisor is not None
    # CAUTION: this function expects polynomials to follow the opposite convention at decoding: the terms must go from the biggest to lowest degree (while most other functions here expect a list from lowest to biggest degree). eg: 1 + 2x + 5x^2 = [5, 2, 1], NOT [1, 2, 5]

    cdef ca.array msg_out = ca.clone(uint64_array_template, dividend.shape[0], False)  # Copy the dividend list and pad with 0 where the ecc bytes will be computed
    cdef uint64_t[::1] msg_out_view = msg_out
    msg_out_view[...] = dividend
    # Cache lengths for faster access inside loops
    cdef uint64_t divisor_len = divisor.shape[0]
    #normalizer = divisor[0]  # precomputing for performance
    cdef uint64_t i, j, coef
    for i in range(dividend.shape[0] - (divisor_len - 1)):
        #msg_out.data.as_ulonglongs[i] /= normalizer  # for general polynomial division (when polynomials are non-monic), the usual way of using synthetic division is to divide the divisor g(x) with its leading coefficient (call it a). In this implementation, this means:we need to compute: coef = msg_out[i] / gen[0]. For more infos, see http://en.wikipedia.org/wiki/Synthetic_division
        coef = msg_out.data.as_ulonglongs[i]  # precaching
        if coef != 0:  # log(0) is undefined, so we need to avoid that case explicitly (and it's also a good optimization). In fact if you remove it, it should still work because gf_mul() will take care of the condition. But it's still a good practice to put the condition here.
            for j in range(1, divisor_len):  # in synthetic division, we always skip the first coefficient of the divisior, because it's only used to normalize the dividend coefficient
                if divisor[j] != 0:  # log(0) is undefined
                    msg_out.data.as_ulonglongs[i + j] ^= gf.gf_mul(divisor[j], coef)  # equivalent to the more mathematically correct (but xoring directly is faster): msg_out[i + j] += -divisor[j] * coef

    # The resulting msg_out contains both the quotient and the remainder, the remainder being the size of the divisor (the remainder has necessarily the same degree as the divisor -- not length but degree == length-1 -- since it's what we couldn't divide from the dividend), so we compute the index where this separation is, and return the quotient and remainder.
    cdef uint64_t separator = -(divisor_len - 1)
    cdef ca.array quotient = ca.clone(uint64_array_template, separator, False)
    cdef uint64_t[::1] quotient_view = quotient
    quotient_view[...] = msg_out_view[:separator]
    cdef ca.array remainder = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>msg_out_view.shape[0] - separator), False)
    cdef uint64_t[::1] remainder_view = remainder
    remainder_view[...] = msg_out_view[separator:]
    return quotient, remainder

cpdef ca.array gf_poly_square(GaloisField gf, const uint64_t[::1] poly): # pragma: no cover
    '''Linear time implementation of polynomial squaring. For details, see paper: "A fast software implementation for arithmetic operations in GF (2n)". De Win, E., Bosselaers, A., Vandenberghe, S., De Gersem, P., & Vandewalle, J. (1996, January). In Advances in Cryptology - Asiacrypt'96 (pp. 65-76). Springer Berlin Heidelberg.'''
    assert poly is not None
    cdef uint64_t length = poly.shape[0]
    cdef ca.array out = ca.clone(uint64_array_template, <Py_ssize_t>(length * 2 - 1), True)
    cdef uint64_t i, p, k
    for i in range(length - 1):
        p = poly[i]
        k = i * 2
        if p != 0:
            #out.data.as_ulonglongs[k] = gf.gf_exp.data.as_ulonglongs[(2 * gf.gf_log.data.as_ulonglongs[p]) % field_charac] # not necessary to modulo (2^r)-1 since gf_exp is duplicated up to 510.
            out.data.as_ulonglongs[k] = gf.gf_exp.data.as_ulonglongs[gf.gf_log.data.as_ulonglongs[p] * 2]
        #else: # not necessary since the output is already initialized to an array of 0
            #out.data.as_ulonglongs[k] = 0
    out.data.as_ulonglongs[length * 2 - 2] = gf.gf_exp.data.as_ulonglongs[gf.gf_log.data.as_ulonglongs[poly[length - 1]] * 2]
    if out.data.as_ulonglongs[0] == 0:
        out.data.as_ulonglongs[0] = poly[1] * 2 - 1
    return out

cpdef uint64_t gf_poly_eval(GaloisField gf, const uint64_t[::1] poly, uint64_t x):
    '''Evaluates a polynomial in GF(2^p) given the value for x. This is based on Horner's scheme for maximum efficiency.'''
    assert poly is not None
    cdef uint64_t y = poly[0]
    cdef uint64_t i
    for i in range(1, <uint64_t>poly.shape[0]):
        y = gf.gf_mul(y, x) ^ poly[i]
    return y


################### REED-SOLOMON ENCODING ###################

cpdef ca.array rs_generator_poly(GaloisField gf, uint64_t nsym, uint64_t fcr, uint64_t generator = 2):
    '''Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)'''
    cdef ca.array g = ca.clone(uint64_array_template, 1, False)
    g.data.as_ulonglongs[0] = 1
    cdef uint64_t i
    cdef ca.array t = ca.clone(uint64_array_template, 2, False)
    t.data.as_ulonglongs[0] = 1
    for i in range(nsym):
        t.data.as_ulonglongs[1] = gf.gf_pow(generator, i + fcr)
        g = gf_poly_mul(gf, g, t)
    return g

cpdef list rs_generator_poly_all(GaloisField gf, uint64_t max_nsym, uint64_t fcr, uint64_t generator = 2):
    '''Generate all irreducible generator polynomials up to max_nsym (usually you can use n, the length of the message+ecc). Very useful to reduce processing time if you want to encode using variable schemes and nsym rates.'''
    #g_all = {}  # old approach using a dict
    #g_all[0] = g_all[1] = pa.array('Q', [1])
    g_all = []  # a list of list is potentially faster than using a dict, since it is pre-allocated and a list has less overhead than a dict
    cdef uint64_t nsym
    for nsym in range(max_nsym):
        g_all.append(rs_generator_poly(gf, nsym, fcr, generator))
    return g_all

cpdef ca.array rs_simple_encode_msg(GaloisField gf, const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = 0, uint64_t generator = 2, const uint64_t[::1] gen = None):
    '''Simple Reed-Solomon encoding (mainly an example for you to understand how it works, because it's slower than the inlined function below)'''
    assert msg_in is not None
    if <uint64_t>msg_in.shape[0] + nsym > gf.field_charac:
        raise ValueError(f'Message is too long ({<uint64_t>msg_in.shape[0] + nsym} when max is {gf.field_charac})')
    if gen is None:
        gen = rs_generator_poly(nsym, fcr, generator)

    # Pad the message, then divide it by the irreducible generator polynomial
    cdef ca.array dividend = ca.clone(uint64_array_template, <Py_ssize_t>(msg_in.shape[0] + gen.shape[0] - 1), True)
    cdef uint64_t[::1] dividend_view = dividend
    dividend_view[:msg_in.shape[0]] = msg_in
    _, remainder = gf_poly_div(gf, dividend_view, gen)
    # The remainder is our RS code! Just append it to our original message to get our full codeword (this represents a polynomial of max 256 terms)
    cdef ca.array msg_out = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>msg_in.shape[0] + <uint64_t>len(remainder)), False)
    cdef uint64_t[::1] msg_out_view = msg_out
    msg_out_view[:msg_in.shape[0]] = msg_in
    msg_out_view[msg_in.shape[0]:] = remainder
    # Return the codeword
    return msg_out

cpdef ca.array rs_encode_msg(GaloisField gf, const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = 0, uint64_t generator = 2, const uint64_t[::1] gen = None):
    '''Reed-Solomon main encoding function, using polynomial division (Extended Synthetic Division, the fastest algorithm available to my knowledge), better explained at http://research.swtch.com/field'''
    assert msg_in is not None
    if <uint64_t>msg_in.shape[0] + nsym > gf.field_charac:
        raise ValueError(f'Message is too long ({<uint64_t>msg_in.shape[0] + nsym} when max is {gf.field_charac})')
    if gen is None:
        gen = rs_generator_poly(nsym, fcr, generator)

    cdef ca.array msg_out = ca.clone(uint64_array_template, <Py_ssize_t>(msg_in.shape[0] + gen.shape[0] - 1), True) # init msg_out with the values inside msg_in and pad with len(gen)-1 bytes (which is the number of ecc symbols).
    cdef uint64_t[::1] msg_out_view = msg_out
    msg_out_view[:msg_in.shape[0]] = msg_in

    # Precompute the logarithm of every items in the generator
    cdef ca.array lgen = ca.clone(uint64_array_template, gen.shape[0], False)
    cdef uint64_t i
    for i in range(<uint64_t>gen.shape[0]):
        lgen.data.as_ulonglongs[i] = gf.gf_log.data.as_ulonglongs[gen[i]]

    # Cache lengths for faster access inside loops
    cdef uint64_t msg_in_len = msg_in.shape[0]
    cdef uint64_t gen_len = gen.shape[0]

    # Extended synthetic division main loop
    # Fastest implementation with PyPy (but the Cython version in creedsolo.pyx is about 2x faster)
    cdef uint64_t coef, j
    for i in range(msg_in_len):
        coef = msg_out.data.as_ulonglongs[i]  # Note that it's msg_out here, not msg_in. Thus, we reuse the updated value at each iteration (this is how Synthetic Division works: instead of storing in a temporary register the intermediate values, we directly commit them to the output).
        # coef = gf_mul(msg_out[i], gf_inverse(gen[0]))  # for general polynomial division (when polynomials are non-monic), the usual way of using synthetic division is to divide the divisor g(x) with its leading coefficient (call it a). In this implementation, this means:we need to compute: coef = msg_out[i] / gen[0]
        if coef != 0:  # log(0) is undefined, so we need to manually check for this case. There's no need to check the divisor here because we know it can't be 0 since we generated it.
            lcoef = gf.gf_log.data.as_ulonglongs[coef]  # precaching

            for j in range(1, gen_len):  # in synthetic division, we always skip the first coefficient of the divisior, because it's only used to normalize the dividend coefficient (which is here useless since the divisor, the generator polynomial, is always monic)
                #if gen[j] != 0:  # log(0) is undefined so we need to check that, but it slow things down in fact and it's useless in our case (reed-solomon encoding) since we know that all coefficients in the generator are not 0
                msg_out.data.as_ulonglongs[i + j] ^= gf.gf_exp.data.as_ulonglongs[lcoef + lgen.data.as_ulonglongs[j]] # optimization, equivalent to gf_mul(gen[j], msg_out[i]) and we just substract it to msg_out[i+j] (but since we are in GF256, it's equivalent to an addition and to an XOR). In other words, this is simply a "multiply-accumulate operation"

    # Recopy the original message bytes (overwrites the part where the quotient was computed)
    msg_out_view[:msg_in_len] = msg_in  # equivalent to c = mprime - b, where mprime is msg_in padded with [0]*nsym
    return msg_out


################### REED-SOLOMON DECODING ###################

cpdef ca.array inverted(const uint64_t[::1] msg):
    '''Implements msg[::-1] explicitly to make the library compatible with MicroPython which does not support stepped slices.'''
    assert msg is not None
    cdef ca.array r = ca.clone(uint64_array_template, msg.shape[0], False)
    cdef uint64_t[::1] r_view = r
    r_view[:] = msg[::-1]
    return r

cpdef ca.array rs_calc_syndromes(GaloisField gf, const uint64_t[::1] msg, uint64_t nsym, uint64_t fcr = 0, uint64_t generator = 2):
    '''Given the received codeword msg and the number of error correcting symbols (nsym), computes the syndromes polynomial.
    Mathematically, it's essentially equivalent to a Fourrier Transform (Chien search being the inverse).
    '''
    assert msg is not None
    # Note the "[0] +" : we add a 0 coefficient for the lowest degree (the constant). This effectively shifts the syndrome, and will shift every computations depending on the syndromes (such as the errors locator polynomial, errors evaluator polynomial, etc. but not the errors positions).
    # This is not necessary as anyway syndromes are defined such as there are only non-zero coefficients (the only 0 is the shift of the constant here) and subsequent computations will/must account for the shift by skipping the first iteration (eg, the often seen range(1, n-k+1)), but you can also avoid prepending the 0 coeff and adapt every subsequent computations to start from 0 instead of 1.
    cdef ca.array r = ca.clone(uint64_array_template, <Py_ssize_t>(nsym + 1), False)
    r.data.as_ulonglongs[0] = 0
    cdef uint64_t i
    for i in range(nsym):
        r.data.as_ulonglongs[i + 1] = gf_poly_eval(gf, msg, gf.gf_pow(generator, i + fcr))
    return r

cpdef ca.array rs_correct_errata(GaloisField gf, const uint64_t[::1] msg_in, const uint64_t[::1] synd, const uint64_t[::1] err_pos, uint64_t fcr = 0, uint64_t generator = 2): # err_pos is a list of the positions of the errors/erasures/errata
    '''Forney algorithm, computes the values (error magnitude) to correct the input message.'''
    assert msg_in is not None and synd is not None and err_pos is not None
    cdef ca.array msg = ca.clone(uint64_array_template, msg_in.shape[0], False)
    cdef uint64_t[::1] msg_view = msg
    msg_view[...] = msg_in
    # calculate errata locator polynomial to correct both errors and erasures (by combining the errors positions given by the error locator polynomial found by BM with the erasures positions given by caller)
    cdef ca.array coef_pos = ca.clone(uint64_array_template, err_pos.shape[0], False)  # need to convert the positions to coefficients degrees for the errata locator algo to work (eg: instead of [0, 1, 2] it will become [len(msg)-1, len(msg)-2, len(msg) -3]), also can use a _bytearray() because length of message - and hence length of coef_pos - will always be constrained by the Galois Field.
    cdef uint64_t i
    for i in range(<uint64_t>err_pos.shape[0]):
        coef_pos.data.as_ulonglongs[i] = <uint64_t>msg_view.shape[0] - 1 - <uint64_t>err_pos[i]
    cdef ca.array err_loc = rs_find_errata_locator(gf, coef_pos, generator)
    # calculate errata evaluator polynomial (often called Omega or Gamma in academic papers)
    cdef ca.array err_eval = inverted(rs_find_error_evaluator(gf, inverted(synd), err_loc, <uint64_t>len(err_loc) - 1))
    # Second part of Chien search to get the error location polynomial x from the error positions in err_pos (the roots of the error locator polynomial, ie, where it evaluates to 0)
    cdef ca.array x = ca.clone(uint64_array_template, len(coef_pos), False)  # will store the position of the errors
    cdef uint64_t l
    for i in range(<uint64_t>len(coef_pos)):
        l = gf.field_charac - coef_pos.data.as_ulonglongs[i]
        x.data.as_ulonglongs[i] = gf.gf_pow(generator, -l)

    # Forney algorithm: compute the magnitudes
    cdef ca.array e = ca.clone(uint64_array_template, len(msg), True)  # will store the values that need to be corrected (substracted) to the message containing errors. This is sometimes called the error magnitude polynomial.
    cdef uint64_t x_len = len(x)
    cdef uint64_t x_i_inv, err_loc_prime, j, y, magnitude
    for i in range(x_len):

        x_i_inv = gf.gf_inverse(x.data.as_ulonglongs[i])

        # Compute the formal derivative of the error locator polynomial (see Blahut, Algebraic codes for data transmission, pp 196-197).
        # the formal derivative of the errata locator is used as the denominator of the Forney Algorithm, which simply says that the ith error value is given by error_evaluator(gf_inverse(Xi)) / error_locator_derivative(gf_inverse(Xi)). See Blahut, Algebraic codes for data transmission, pp 196-197.
        #err_loc_prime_tmp = pa.array('Q')
        #for j in range(x_len):
        #    if j != i:
        #        err_loc_prime_tmp.append(gf_sub(1, gf.gf_mul(x_i_inv, x[j])))
        # compute the product, which is the denominator of the Forney algorithm (errata locator derivative)
        #err_loc_prime = 1
        #for coef in err_loc_prime_tmp:
        #    err_loc_prime = gf_mul(err_loc_prime, coef)
        # equivalent to: err_loc_prime = functools.reduce(gf.gf_mul, err_loc_prime_tmp, 1)

        # Alternative but faster way to compute the formal derivative of the error locator polynomial
        err_loc_prime = 1
        for j in range(x_len):
            if j != i:
                err_loc_prime = gf.gf_mul(err_loc_prime, gf_sub(1, gf.gf_mul(x_i_inv, x.data.as_ulonglongs[j])))

        # Test if we could find the errata locator, else we raise an Exception (because else since we divide y by err_loc_prime to compute the magnitude, we will get a ZeroDivisionError exception otherwise)
        if err_loc_prime == 0:
            raise ReedSolomonError('Decoding failed: Forney algorithm could not properly detect where the errors are located (errata locator prime is 0).')

        # Compute y (evaluation of the errata evaluator polynomial)
        # This is a more faithful translation of the theoretical equation contrary to the old forney method. Here it is exactly copy/pasted from the included presentation decoding_rs.pdf: Yl = omega(Xl.inverse()) / prod(1 - Xj*Xl.inverse()) for j in len(X) (in the paper it's for j in s, but it's useless when len(X) < s because we compute neutral terms 1 for nothing, and wrong when correcting more than s erasures or erasures+errors since it prevents computing all required terms).
        # Thus here this method works with erasures too because firstly we fixed the equation to be like the theoretical one (don't know why it was modified in _old_forney(), if it's an optimization, it doesn't enhance anything), and secondly because we removed the product bound on s, which prevented computing errors and erasures above the s=(n-k)//2 bound.
        y = gf_poly_eval(gf, inverted(err_eval), x_i_inv)  # numerator of the Forney algorithm (errata evaluator evaluated)
        y = gf.gf_mul(gf.gf_pow(x.data.as_longlongs[i], 1 - <int64_t>fcr), y)  # adjust to fcr parameter

        # Compute the magnitude
        magnitude = gf.gf_div(y, err_loc_prime)  # magnitude value of the error, calculated by the Forney algorithm (an equation in fact): dividing the errata evaluator with the errata locator derivative gives us the errata magnitude (ie, value to repair) the ith symbol
        e.data.as_ulonglongs[err_pos[i]] = magnitude  # store the magnitude for this error into the magnitude polynomial

    # Apply the correction of values to get our message corrected! (note that the ecc bytes also gets corrected!)
    # (this isn't the Forney algorithm, we just apply the result of decoding here)
    msg = gf_poly_add(msg, e)  # equivalent to Ci = Ri - Ei where Ci is the correct message, Ri the received (senseword) message, and Ei the errata magnitudes (minus is replaced by XOR since it's equivalent in GF(2^p)). So in fact here we substract from the received message the errors magnitude, which logically corrects the value to what it should be.
    return msg

cpdef ca.array rs_find_error_locator(GaloisField gf, const uint64_t[::1] synd, uint64_t nsym, const uint64_t[::1] erase_loc = None, uint64_t erase_count = 0):
    '''Find error/errata locator and evaluator polynomials with Berlekamp-Massey algorithm'''
    assert synd is not None
    # The idea is that BM will iteratively estimate the error locator polynomial.
    # To do this, it will compute a Discrepancy term called Delta, which will tell us if the error locator polynomial needs an update or not
    # (hence why it's called discrepancy: it tells us when we are getting off board from the correct value).

    # Init the polynomials
    cdef ca.array err_loc, old_loc
    cdef uint64_t[::1] err_loc_view, old_loc_view
    cdef uint64_t i
    if erase_loc is not None:  # if the erasure locator polynomial is supplied, we init with its value, so that we include erasures in the final locator polynomial
        err_loc = ca.clone(uint64_array_template, erase_loc.shape[0], False)
        err_loc_view = err_loc
        err_loc_view[:] = err_loc
        old_loc = ca.clone(uint64_array_template, erase_loc.shape[0], False)
        old_loc_view = old_loc
        old_loc_view[:] = erase_loc
    else:
        err_loc = ca.clone(uint64_array_template, 1, False)  # This is the main variable we want to fill, also called Sigma in other notations or more formally the errors/errata locator polynomial.
        err_loc.data.as_ulonglongs[0] = 1
        old_loc = ca.clone(uint64_array_template, 1, False)  # BM is an iterative algorithm, and we need the errata locator polynomial of the previous iteration in order to update other necessary variables.
        old_loc.data.as_ulonglongs[0] = 1
    #L = 0 # update flag variable, not needed here because we use an alternative equivalent way of checking if update is needed (but using the flag could potentially be faster depending on if using length(list) is taking linear time in your language, here in Python it's constant so it's as fast.

    # Fix the syndrome shifting: when computing the syndrome, some implementations may prepend a 0 coefficient for the lowest degree term (the constant). This is a case of syndrome shifting, thus the syndrome will be bigger than the number of ecc symbols (I don't know what purpose serves this shifting). If that's the case, then we need to account for the syndrome shifting when we use the syndrome such as inside BM, by skipping those prepended coefficients.
    # Another way to detect the shifting is to detect the 0 coefficients: by definition, a syndrome does not contain any 0 coefficient (except if there are no errors/erasures, in this case they are all 0). This however doesn't work with the modified Forney syndrome, which set to 0 the coefficients corresponding to erasures, leaving only the coefficients corresponding to errors.
    cdef uint64_t synd_shift = 0
    if <uint64_t>synd.shape[0] > nsym:
        synd_shift = synd.shape[0] - nsym

    cdef uint64_t k, delta, j, err_loc_len
    cdef ca.array new_loc
    for i in range(nsym - erase_count):  # generally: nsym-erase_count == len(synd), except when you input a partial erase_loc and using the full syndrome instead of the Forney syndrome, in which case nsym-erase_count is more correct (len(synd) will fail badly with IndexError).
        if erase_loc is not None: # if an erasures locator polynomial was provided to init the errors locator polynomial, then we must skip the FIRST erase_count iterations (not the last iterations, this is very important!)
            k = erase_count + i + synd_shift
        else:  # if erasures locator is not provided, then either there's no erasures to account or we use the Forney syndromes, so we don't need to use erase_count nor erase_loc (the erasures have been trimmed out of the Forney syndromes).
            k = i + synd_shift

        # Compute the discrepancy Delta
        # Here is the close-to-the-books operation to compute the discrepancy Delta: it's a simple polynomial multiplication of error locator with the syndromes, and then we get the Kth element.
        #delta = gf_poly_mul(inverted(err_loc), synd)[K] # theoretically it should be inverted(gf_poly_add(inverted(synd), [1])) instead of just synd, but it seems it's not absolutely necessary to correctly decode.
        # But this can be optimized: since we only need the Kth element, we don't need to compute the polynomial multiplication for any other element but the Kth. Thus to optimize, we compute the polymul only at the item we need, skipping the rest (avoiding a nested loop, thus we are linear time instead of quadratic).
        # This optimization is actually described in several figures of the book "Algebraic codes for data transmission", Blahut, Richard E., 2003, Cambridge university press.
        delta = synd[k]
        err_loc_len = len(err_loc)
        for j in range(1, <uint64_t>len(err_loc)):  # range 1:256 is important: if you use range 0:255, if the last byte of the ecc symbols is corrupted, it won't be correctable! You need to use the range 1,256 to include this last byte.
            delta ^= gf.gf_mul(err_loc.data.as_ulonglongs[err_loc_len - (j + 1)], synd[k - j]) # delta is also called discrepancy. Here we do a partial polynomial multiplication (ie, we compute the polynomial multiplication only for the term of degree K). Should be equivalent to brownanrs.polynomial.mul_at().
        #print('delta', k, delta, list(gf_poly_mul(inverted(err_loc), synd))) # debugline

        # Shift polynomials to compute the next degree
        ca.resize(old_loc, <Py_ssize_t>len(old_loc) + 1)
        old_loc.data.as_longlongs[len(old_loc) - 1] = 0

        # Iteratively estimate the errata locator and evaluator polynomials
        if delta != 0:  # Update only if there's a discrepancy
            if len(old_loc) > len(err_loc):  # Rule B (rule A is implicitly defined because rule A just says that we skip any modification for this iteration)
            #if 2 * l <= k + erase_count:  # equivalent to len(old_loc) > len(err_loc), as long as L is correctly computed
                # Computing errata locator polynomial Sigma
                new_loc = gf_poly_scale(gf, old_loc, delta)
                old_loc = gf_poly_scale(gf, err_loc, gf.gf_inverse(delta))  # effectively we are doing err_loc * 1/delta = err_loc // delta
                err_loc = new_loc
                # Update the update flag
                #l = k - l  # the update flag L is tricky: in Blahut's schema, it's mandatory to use `L = K - L - erase_count` (and indeed in a previous draft of this function, if you forgot to do `- erase_count` it would lead to correcting only 2*(errors+erasures) <= (n-k) instead of 2*errors+erasures <= (n-k)), but in this latest draft, this will lead to a wrong decoding in some cases where it should correctly decode! Thus you should try with and without `- erase_count` to update L on your own implementation and see which one works OK without producing wrong decoding failures.

            # Update with the discrepancy
            err_loc = gf_poly_add(err_loc, gf_poly_scale(gf, old_loc, delta))

    # Check if the result is correct, that there's not too many errors to correct
    cdef uint64_t[::1] new_loc_view
    for i in range(<uint64_t>len(err_loc)):  # drop leading 0s, else errs will not be of the correct size. This does not use functional closures (ie, lambdas), which Cython and JIT compilers cannot optimize, such as `err_loc = list(itertools.dropwhile(lambda x: x == 0, err_loc))`
        if err_loc.data.as_ulonglongs[i] != 0:
            new_loc = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>len(err_loc) - i), False)
            new_loc_view = new_loc
            err_loc_view = err_loc
            new_loc_view[:] = err_loc_view[i:]
            err_loc = new_loc
            break

    cdef uint64_t errs = <uint64_t>len(err_loc) - 1  # -1 because range is 1:256, it's offset by 1, it's not 0:255, hence the length would be overestimated without -1
    if errs * 2 > nsym + erase_count:  # failure if we have too many erratas for the Singleton Bound.
        raise ReedSolomonError('Too many errors to correct')

    # # Return result
    return err_loc

cpdef ca.array rs_find_errata_locator(GaloisField gf, const uint64_t[::1] e_pos, uint64_t generator = 2):
    '''Compute the erasures/errors/errata locator polynomial from the erasures/errors/errata positions (the positions must be relative to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to "h_ll_ worldxxxxxxxxx" with xxxxxxxxx being the ecc of length n-k=9, here the string positions are [1, 4], but the coefficients are reversed since the ecc characters are placed as the first coefficients of the polynomial, thus the coefficients of the erased characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified as an argument.'''
    assert e_pos is not None
    # See: http://ocw.usu.edu/Electrical_and_Computer_Engineering/Error_Control_Coding/lecture7.pdf and Blahut, Richard E. "Transform techniques for error control codes." IBM Journal of Research and development 23.3 (1979): 299-315. http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.600&rep=rep1&type=pdf and also a MatLab implementation here: http://www.mathworks.com/matlabcentral/fileexchange/23567-reed-solomon-errors-and-erasures-decoder/content//RS_E_E_DEC.m
    cdef ca.array e_loc = ca.clone(uint64_array_template, 1, False) # just to init because we will multiply, so it must be 1 so that the multiplication starts correctly without nulling any term
    e_loc.data.as_ulonglongs[0] = 1
    # erasures_loc is very simple to compute: erasures_loc = prod(1 - x * alpha ** i) for i in erasures_pos and where alpha is the alpha chosen to evaluate polynomials (here in this library it's gf(3)). To generate c*x where c is a constant, we simply generate a Polynomial([c, 0]) where 0 is the constant and c is positionned to be the coefficient for x^1.
    cdef uint64_t i
    cdef ca.array one = ca.clone(uint64_array_template, 1, False)
    one.data.as_ulonglongs[0] = 1
    cdef ca.array t = ca.clone(uint64_array_template, 2, False)
    t.data.as_ulonglongs[1] = 0
    for i in range(<uint64_t>e_pos.shape[0]):
        t.data.as_ulonglongs[0] = gf.gf_pow(generator, e_pos[i])
        e_loc = gf_poly_mul(gf, e_loc, gf_poly_add(one, t))
    return e_loc

cpdef ca.array rs_find_error_evaluator(GaloisField gf, const uint64_t[::1] synd, const uint64_t[::1] err_loc, uint64_t nsym):
    '''Compute the error (or erasures if you supply sigma=erasures locator polynomial, or errata) evaluator polynomial Omega from the syndrome and the error/erasures/errata locator Sigma. Omega is already computed at the same time as Sigma inside the Berlekamp-Massey implemented above, but in case you modify Sigma, you can recompute Omega afterwards using this method, or just ensure that Omega computed by BM is correct given Sigma.'''
    assert synd is not None and err_loc is not None
    # Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)
    #_, remainder = gf_poly_div(gf_poly_mul(synd, err_loc), pa.array('Q', [1] + [0] * (nsym + 1)))  # first multiply syndromes * errata_locator, then do a polynomial division to truncate the polynomial to the required length

    # Faster way that is equivalent
    cdef ca.array remainder = gf_poly_mul(gf, synd, err_loc)  # first multiply the syndromes with the errata locator polynomial
    cdef ca.array r = ca.clone(uint64_array_template, <Py_ssize_t>(nsym + 1), False)
    cdef uint64_t[::1] r_view = r
    cdef uint64_t[::1] remainder_view = remainder
    r_view[:] = remainder_view[<uint64_t>remainder_view.shape[0] - (nsym + 1):]
    return r

cpdef ca.array rs_find_errors(GaloisField gf, const uint64_t[::1] err_loc, uint64_t nmess, uint64_t generator = 2):
    '''Find the roots (ie, where evaluation = zero) of error polynomial by smart bruteforce trial. This is a faster form of chien search, processing only useful coefficients (the ones in the messages) instead of the whole 2^8 range. Besides the speed boost, this also allows to fix a number of issue: correctly decoding when the last ecc byte is corrupted, and accepting messages of length n > 2^8.'''
    assert err_loc is not None
    # nmess = length of whole codeword (message + ecc symbols)
    cdef ca.array err_pos = ca.clone(uint64_array_template, nmess, False)
    cdef uint64_t err_pos_size = 0
    cdef uint64_t i
    for i in range(nmess):  # normally we should try all 2^8 possible values, but here we optimize to just check the interesting symbols
        if gf_poly_eval(gf, err_loc, gf.gf_pow(generator, nmess - 1 - i)) == 0:  # It's a 0? Bingo, it's a root of the error locator polynomial, in other terms this is the location of an error
            err_pos.data.as_ulonglongs[err_pos_size] = i
            err_pos_size += 1
    ca.resize(err_pos, err_pos_size)
    # Sanity check: the number of errors/errata positions found should be exactly the same as the length of the errata locator polynomial
    cdef uint64_t errs = <uint64_t>err_loc.shape[0] - 1  # compute the exact number of errors/errata that this error locator should find
    if <uint64_t>len(err_pos) != errs:
        # TODO: to decode messages+ecc with length n > 255, we may try to use a bruteforce approach: the correct positions ARE in the final array j, but the problem is because we are above the Galois Field's range, there is a wraparound so that for example if j should be [0, 1, 2, 3], we will also get [255, 256, 257, 258] (because 258 % 255 == 3, same for the other values), so we can't discriminate. The issue is that fixing any errs_nb errors among those will always give a correct output message (in the sense that the syndrome will be all 0), so we may not even be able to check if that's correct or not, so I'm not sure the bruteforce approach may even be possible.
        raise ReedSolomonError('Too many (or few) errors found by Chien Search for the errata locator polynomial!')
    return err_pos

cpdef ca.array rs_forney_syndromes(GaloisField gf, const uint64_t[::1] synd, const uint64_t[::1] pos, uint64_t nmess, uint64_t generator = 2):
    # Compute Forney syndromes, which computes a modified syndromes to compute only errors (erasures are trimmed out). Do not confuse this with Forney algorithm, which allows to correct the message based on the location of errors.
    assert synd is not None and pos is not None
    cdef ca.array erase_pos_reversed = ca.clone(uint64_array_template, pos.shape[0], False)  # prepare the coefficient degree positions (instead of the erasures positions)
    cdef uint64_t i
    for i in range(<uint64_t>pos.shape[0]):
        erase_pos_reversed.data.as_ulonglongs[i] = nmess - 1 - <uint64_t>pos[i]
    # Optimized method, all operations are inlined
    cdef ca.array fsynd = ca.clone(uint64_array_template, <Py_ssize_t>(synd.shape[0] - 1), False)  # make a copy and trim the first coefficient which is always 0 by definition
    cdef uint64_t[::1] fsynd_view = fsynd
    fsynd_view[...] = synd[1:]
    cdef uint64_t x, j
    for i in range(<uint64_t>pos.shape[0]):
        pass
        x = gf.gf_pow(generator, erase_pos_reversed.data.as_ulonglongs[i])
        for j in range(<uint64_t>len(fsynd) - 1):
            fsynd.data.as_ulonglongs[j] = gf.gf_mul(fsynd.data.as_ulonglongs[j], x) ^ fsynd.data.as_ulonglongs[j + 1]
        #fsynd.pop()  # useless? it doesn't change the results of computations to leave it there

    # Theoretical way of computing the modified Forney syndromes: fsynd = (erase_loc * synd) % x^(n-k) -- although the trimming by using x^(n-k) is maybe not necessary as many books do not even mention it (and it works without trimming)
    # See Shao, H. M., Truong, T. K., Deutsch, L. J., & Reed, I. S. (1986, April). A single chip VLSI Reed-Solomon decoder. In Acoustics, Speech, and Signal Processing, IEEE International Conference on ICASSP'86. (Vol. 11, pp. 2151-2154). IEEE.ISO 690
    #erase_loc = rs_find_errata_locator(gf, erase_pos_reversed, generator)  # computing the erasures locator polynomial
    #fsynd = gf_poly_mul(gf, inverted(erase_loc), synd[1:])  # then multiply with the syndrome to get the untrimmed forney syndrome
    #fsynd = fsynd[len(pos):]  # then trim the first erase_pos coefficients which are useless. Seems to be not necessary, but this reduces the computation time later in BM (thus it's an optimization).

    return fsynd

cpdef tuple rs_correct_msg(GaloisField gf, const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = 0, uint64_t generator = 2, const uint64_t[::1] erase_pos = None, bool only_erasures = False):
    '''Reed-Solomon main decoding function'''
    assert msg_in is not None
    if <uint64_t>msg_in.shape[0] > gf.field_charac:
        # Note that it is in fact possible to encode/decode messages that are longer than field_charac, but because this will be above the field, this will generate more error positions during Chien Search than it should, because this will generate duplicate values, which should normally be prevented thank's to the prime polynomial reduction (eg, because it can't discriminate between error at position 1 or 256, both being exactly equal under galois field 2^8). So it's really not advised to do it, but it's possible (but then you're not guaranted to be able to correct any error/erasure on symbols with a position above the length of field_charac -- if you really need a bigger message without chunking, then you should better enlarge c_exp so that you get a bigger field).
        raise ValueError(f'Message is too long ({msg_in.shape[0]} when max is {gf.field_charac})')

    cdef ca.array msg_out = ca.clone(uint64_array_template, msg_in.shape[0], False)  # copy of message
    cdef uint64_t[::1] msg_out_view = msg_out
    msg_out_view[:] = msg_in
    # erasures: set them to null bytes for easier decoding (but this is not necessary, they will be corrected anyway, but debugging will be easier with null bytes because the error locator polynomial values will only depend on the errors locations, not their values)
    cdef uint64_t i
    if erase_pos is None:
        erase_pos = ca.clone(uint64_array_template, 0, False)
    else:
        for i in range(<uint64_t>erase_pos.shape[0]):
            msg_out.data.as_ulonglongs[erase_pos[i]] = 0
    # check if there are too many erasures to correct (beyond the Singleton bound)
    if <uint64_t>len(erase_pos) > nsym:
        raise ReedSolomonError('Too many erasures to correct')
    # prepare the syndrome polynomial using only errors (ie: errors = characters that were either replaced by null byte or changed to another character, but we don't know their positions)
    cdef ca.array synd = rs_calc_syndromes(gf, msg_out, nsym, fcr, generator)
    # check if there's any error/erasure in the input codeword. If not (all syndromes coefficients are 0), then just return the codeword as-is.
    cdef bool flag = True
    for i in range(<uint64_t>len(synd)):
        if synd.data.as_ulonglongs[i]:
            flag = False
            break
    cdef ca.array content, ecc, wrong_pos
    cdef uint64_t[::1] content_view, ecc_view, wrong_pos_view
    if flag:
        content = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>msg_out_view.shape[0] - nsym), False)
        content_view = content
        content_view[:] = msg_out_view[:-nsym]
        ecc = ca.clone(uint64_array_template, nsym, False)
        ecc_view = ecc
        ecc_view[:] = msg_out_view[-nsym:]
        wrong_pos = ca.clone(uint64_array_template, erase_pos.shape[0], False)
        wrong_pos_view = wrong_pos
        wrong_pos_view[:] = erase_pos
        return content, ecc, wrong_pos # no errors

    # Find errors locations
    cdef ca.array err_pos, fsynd, err_loc
    if only_erasures:
        err_pos = ca.clone(uint64_array_template, 0, False)
    else:
        # compute the Forney syndromes, which hide the erasures from the original syndrome (so that BM will just have to deal with errors, not erasures)
        fsynd = rs_forney_syndromes(gf, synd, erase_pos, len(msg_out), generator)
        # compute the error locator polynomial using Berlekamp-Massey
        err_loc = rs_find_error_locator(gf, fsynd, nsym, erase_count = erase_pos.shape[0])
        # locate the message errors using Chien search (or bruteforce search)
        err_pos = rs_find_errors(gf, inverted(err_loc), len(msg_out), generator)
        if err_pos is None:
            raise ReedSolomonError('Could not locate error')

    # Find errors values and apply them to correct the message
    # compute errata evaluator and errata magnitude polynomials, then correct errors and erasures
    cdef uint64_t[::1] err_pos_view = err_pos
    wrong_pos = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>erase_pos.shape[0] + <uint64_t>len(err_pos)), False)
    wrong_pos_view = wrong_pos
    wrong_pos_view[:erase_pos.shape[0]] = erase_pos
    wrong_pos_view[erase_pos.shape[0]:] = err_pos_view
    msg_out = rs_correct_errata(gf, msg_out, synd, wrong_pos, fcr, generator)  # note that we here use the original syndrome, not the forney syndrome (because we will correct both errors and erasures, so we need the full syndrome)
    # check if the final message is fully repaired
    synd = rs_calc_syndromes(gf, msg_out, nsym, fcr, generator)
    flag = True
    for i in range(<uint64_t>len(synd)):
        if synd.data.as_ulonglongs[i] != 0:
            flag = False
            break
    if not flag:
        raise ReedSolomonError('Could not correct message')
    # return the successfully decoded message
    msg_out_view = msg_out
    content = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>msg_out_view.shape[0] - nsym), False)
    content_view = content
    content_view[:] = msg_out_view[:-nsym]
    ecc = ca.clone(uint64_array_template, nsym, False)
    ecc_view = ecc
    ecc_view[:] = msg_out_view[-nsym:]
    return content, ecc, wrong_pos  # also return the corrected ecc block so that the user can check(), and the position of errors to allow for adaptive bitrate algorithm to check how the number of errors vary

# cpdef tuple rs_correct_msg_nofsynd(const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = 0, uint64_t generator = 2, const uint64_t[::1] erase_pos = None, bool only_erasures = False):
#     '''Reed-Solomon main decoding function, without using the modified Forney syndromes'''
#     assert msg_in is not None
#     if <uint64_t>msg_in.shape[0] > field_charac:
#         raise ValueError(f'Message is too long ({msg_in.shape[0]} when max is {field_charac})')

#     cdef ca.array msg_out = ca.clone(uint64_array_template, msg_in.shape[0], False)  # copy of message
#     cdef uint64_t[::1] msg_out_view = msg_out
#     msg_out_view[:] = msg_in
#     # erasures: set them to null bytes for easier decoding (but this is not necessary, they will be corrected anyway, but debugging will be easier with null bytes because the error locator polynomial values will only depend on the errors locations, not their values)
#     cdef uint64_t i
#     if erase_pos is None:
#         erase_pos = ca.clone(uint64_array_template, 0, False)
#     else:
#         for i in range(<uint64_t>erase_pos.shape[0]):
#             msg_out.data.as_ulonglongs[erase_pos[i]] = 0
#     # check if there are too many erasures
#     if <uint64_t>erase_pos.shape[0] > nsym:
#         raise ReedSolomonError(Too many erasures to correct')
#     # prepare the syndrome polynomial using only errors (ie: errors = characters that were either replaced by null byte or changed to another character, but we don't know their positions)
#     cdef ca.array synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
#     # check if there's any error/erasure in the input codeword. If not (all syndromes coefficients are 0), then just return the codeword as-is.
#     cdef bool flag = True
#     for i in range(<uint64_t>len(synd)):
#         if synd.data.as_ulonglongs[i]:
#             flag = False
#             break
#     cdef ca.array content, ecc, wrong_pos
#     cdef uint64_t content_view, ecc_view, wrong_pos_view
#     if flag:
#         content = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>msg_out_view.shape[0] - nsym), False)
#         content_view = content
#         content_view[:] = msg_out_view[:-nsym]
#         ecc = ca.clone(uint64_array_template, nsym, False)
#         ecc_view = ecc
#         ecc_view[:] = msg_out_view[-nsym:]
#         wrong_pos = ca.clone(uint64_array_template, erase_pos.shape[0], False)
#         wrong_pos_view = wrong_pos
#         wrong_pos_view[:] = erase_pos
#         return content, ecc, wrong_pos # no errors

#     # prepare erasures locator and evaluator polynomials
#     cdef uint64_t erase_count = erase_pos.shape[0]
#     cdef uint64_t msg_out_len = len(msg_out)  # cache to avoid recalculations inside loop
#     cdef ca.array erase_pos_reversed = ca.clone(uint64_array_template, erase_pos.shape[0], False)
#     for i in range(<uint64_t>erase_pos.shape[0]):
#         erase_pos_reversed.data.as_ulonglongs[i] = msg_out_len - 1 - <uint64_t>erase_pos[i]
#     cdef ca.array erase_loc = rs_find_errata_locator(erase_pos_reversed, generator)
#     #cdef ca.array erase_eval = rs_find_error_evaluator(inverted(synd), erase_loc, <uint64_t>len(erase_loc) - 1)

#     # prepare errors/errata locator polynomial
#     cdef ca.array err_loc
#     if only_erasures:
#         err_loc = inverted(erase_loc)
#         #err_eval = inverted(erase_eval)
#     else:
#         err_loc = rs_find_error_locator(synd, nsym, erase_loc, erase_count)
#         err_loc = inverted(err_loc)
#         #err_eval = inverted(rs_find_error_evaluator(inverted(synd), inverted(err_loc), <uint64_t>len(err_loc) - 1))  # find error/errata evaluator polynomial (not really necessary since we already compute it at the same time as the error locator poly in BM)

#     # locate the message errors
#     cdef ca.array err_pos = rs_find_errors(err_loc, len(msg_out), generator)  # find the roots of the errata locator polynomial (ie: the positions of the errors/errata)
#     if err_pos is None:
#         raise ReedSolomonError('Could not locate error')

#     # compute errata evaluator and errata magnitude polynomials, then correct errors and erasures
#     msg_out = rs_correct_errata(msg_out, synd, err_pos, fcr, generator)
#     # check if the final message is fully repaired
#     synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
#     flag = True
#     for i in range(<uint64_t>len(synd)):
#         if synd.data.as_ulonglongs[i]:
#             flag = False
#             break
#     if not flag:
#         raise ReedSolomonError('Could not correct message')
#     # return the successfully decoded message
#     msg_out_view = msg_out
#     content = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>msg_out_view.shape[0] - nsym), False)
#     content_view = content
#     content_view[:] = msg_out_view[:-nsym]
#     ecc = ca.clone(uint64_array_template, nsym, False)
#     ecc_view = ecc
#     ecc_view[:] = msg_out_view[-nsym:]
#     cdef uint64_t[::1] err_pos_view = err_pos
#     wrong_pos = ca.clone(uint64_array_template, <Py_ssize_t>(<uint64_t>erase_pos.shape[0] + <uint64_t>len(err_pos)), False)
#     wrong_pos_view = wrong_pos
#     wrong_pos_view[:erase_pos.shape[0]] = erase_pos
#     wrong_pos_view[erase_pos.shape[0]:] = err_pos_view
#     return content, ecc, wrong_pos  # also return the corrected ecc block so that the user can check(), and the position of errors to allow for adaptive bitrate algorithm to check how the number of errors vary

cpdef bool rs_check(GaloisField gf, const uint64_t[::1] msg, uint64_t nsym, uint64_t fcr = 0, uint64_t generator = 2):
    '''Returns true if the message + ecc has no error of false otherwise (may not always catch a wrong decoding or a wrong message, particularly if there are too many errors -- above the Singleton bound --, but it usually does)'''
    assert msg is not None
    cdef ca.array synd = rs_calc_syndromes(gf, msg, nsym, fcr, generator)
    for i in range(<uint64_t>len(synd)):
        if synd.data.as_ulonglongs[i]:
            return False
    return True


#===================================================================================================
# API
#===================================================================================================
cdef class RSCodec:
    '''
    A Reed Solomon encoder/decoder. After initializing the object, use ``encode`` to encode a
    (byte)string to include the RS correction code, and pass such an encoded (byte)string to
    ``decode`` to extract the original message (if the number of errors allows for correct decoding).
    The ``nsym`` argument is the length of the correction code, and it determines the number of
    error bytes (if I understand this correctly, half of ``nsym`` is correctable)
    '''
    '''
    Modifications by rotorgit 2/3/2015:
    Added support for US FAA ADSB UAT RS FEC, by allowing user to specify
    different primitive polynomial and non-zero first consecutive root (fcr).
    For UAT/ADSB use, set fcr=120 and prim=0x187 when instantiating
    the class; leaving them out will default for previous values (0 and
    0x11d)
    '''

    def __init__(self, nsym: int = 10, nsize: int = 255, fcr: int = 0, prim: int = 0x11d, generator: int = 2, c_exp=8, single_gen: bool = True) -> None:
        '''Initialize the Reed-Solomon codec. Note that different parameters change the internal values (the ecc symbols, look-up table values, etc) but not the output result (whether your message can be repaired or not, there is no influence of the parameters).
        nsym : number of ecc symbols (you can repair nsym/2 errors and nsym erasures.
        nsize : maximum length of each chunk. If higher than 255, will use a higher Galois Field, but the algorithm's complexity and computational cost will raise quadratically...
        single_gen : if you want to use the same RSCodec for different nsym parameters (but nsize the same), then set single_gen=False. This is only required for encoding with various number of ecc symbols, as for decoding this is always possible even if single_gen=True.
        '''

        # Auto-setup if galois field or message length is different than default (exponent 8)
        cdef uint64_t n
        if nsize > 255 and c_exp <= 8:  # nsize (chunksize) is larger than the galois field, we resize the galois field
            # Get the next closest power of two
            n = nsize
            c_exp = 0
            while n:
                n >>= 1
                c_exp += 1
        if c_exp != 8 and prim == 0x11d:  # prim was not correctly defined, find one
            prim = <int>find_prime_polys(generator, c_exp, True, True).data.as_ulonglongs[0]
            if nsize == 255:  # resize chunk size if not set
                nsize = (1 << c_exp) - 1
        if nsym >= nsize:
            raise ValueError('ECC symbols must be strictly less than the total message length (nsym < nsize).')

        # Memorize variables
        self.nsym = nsym  # number of ecc symbols (ie, the repairing rate will be r=(nsym/2)/nsize, so for example if you have nsym=5 and nsize=10, you have a rate r=0.25, so you can correct up to 0.25% errors (or exactly 2 symbols out of 10), and 0.5% erasures (5 symbols out of 10).
        self.nsize = nsize  # maximum length of one chunk (ie, message + ecc symbols after encoding, for the message alone it's nsize-nsym)
        self.fcr = fcr  # first consecutive root, can be any value between 0 and (2**c_exp)-1
        self.prim = prim  # prime irreducible polynomial, use find_prime_polys() to find a prime poly
        self.generator = generator  # generator integer, must be prime
        self.c_exp = c_exp  # exponent of the field's characteristic. This both defines the maximum value per symbol and the maximum length of one chunk. By default it's GF(2^8), do not change if you're not sure what it means.
        self.single_gen = single_gen

        # Initialize the look-up tables for easy and quick multiplication/division
        self.gf = GaloisField(prim, generator, c_exp)
        # Precompute the generator polynomials
        if single_gen:
            self.gen = rs_generator_poly(self.gf, nsym, fcr, generator)
        else:
            self.gens = rs_generator_poly_all(self.gf, nsize, fcr, generator)

    def chunk(self, data: Sequence[int], chunk_size: int) -> Generator[Sequence[int], None, None]:
        '''Split a long message into chunks
        DEPRECATED: inlined alternate form so that we can preallocate arrays and hence get faster results with JIT compilers such as PyPy.'''
        cdef uint64_t i
        for i in range(len(data), <uint64_t>chunk_size):
            # alternative chunking form: j = i * chunk_size; chunk = data[j:j + chunk_size]
            # Split the long message in a chunk
            chunk = data[i:i + chunk_size]
            yield chunk

    def encode(self, data: Iterable[int], nsym: Optional[int] = None) -> Annotated[pa.array[int], 'Q']:
        '''Encode a message (ie, add the ecc symbols) using Reed-Solomon, whatever the length of the message because we use chunking
        Optionally, can set nsym to encode with a different number of error correction symbols, but RSCodec must be initialized with single_gen=False first.
        slice_assign=True allows to speed up the loop quite significantly in JIT compilers such as PyPy by preallocating the output bytearray and slice assigning into it, instead of constantly extending an empty bytearray, but this only works in Python 3, not Python 2, hence is disabled by default for retrocompatibility.
        '''

        cdef uint64_t nsize = self.nsize
        cdef uint64_t fcr = self.fcr
        cdef uint64_t generator = self.generator

        cdef uint64_t nsym_c
        if nsym is None:
            nsym_c = self.nsym
        else:
            nsym_c = nsym

        cdef ca.array gen
        if self.single_gen:
            if nsym_c != self.nsym:
                gen = rs_generator_poly(self.gf, nsym_c, fcr, generator)
            else:
                gen = self.gen
        else:
            gen = self.gens[nsym_c]

        cdef ca.array data_array = <ca.array>pa.array('Q', data)
        cdef uint64_t[::1] data_array_view = data_array

        # Calculate chunk size and total number of chunks for looping
        cdef uint64_t chunk_size = nsize - nsym_c
        cdef uint64_t total_chunks = (<uint64_t>len(data_array) - 1) // chunk_size + 1  # need to convert to floats first to get an accurate floating division, or else we assume implicit conversion and it will cause an error on Python 2

        # Preallocate output array
        cdef ca.array enc = ca.clone(uint64_array_template, <Py_ssize_t>(nsize * total_chunks), False)  # pre-allocate array and we will overwrite data in it, much faster than extending
        cdef uint64_t[::1] enc_view = enc
        # Chunking loop
        cdef uint64_t i
        cdef ca.array r
        cdef uint64_t[::1] r_view
        for i in range(total_chunks):
            # Encode this chunk and update a slice of the output bytearray, much more efficient than extending an array constantly
            r = rs_encode_msg(self.gf, data_array_view[chunk_size * i:chunk_size * (i + 1)], nsym_c, fcr, generator, gen)
            r_view = r
            enc_view[nsize * i:nsize * i + r_view.shape[0]] = r_view
            if i == total_chunks - 1:
                ca.resize(enc, <Py_ssize_t>(nsize * i + <uint64_t>r_view.shape[0]))

        return enc

    def decode(self, data: Iterable[int], nsym: Optional[int] = None, erase_pos: Optional[Iterable[int]] = None, only_erasures: bool = False) -> tuple[Annotated[pa.array[int], 'Q'], Annotated[pa.array[int], 'Q'], Annotated[pa.array[int], 'Q']]:
        '''Repair a message, whatever its size is, by using chunking. May return a wrong result if number of errors > nsym because then too many errors to be corrected.
        Note that it returns a couple of vars: the repaired messages, and the repaired messages+ecc (useful for checking).
        Usage: rmes, rmesecc = RSCodec.decode(data).
        Optionally: can specify nsym to decode messages of different parameters, erase_pos with a list of erasures positions to double the number of erasures that can be corrected compared to unlocalized errors, only_erasures boolean to specify if we should only look for erasures, which speeds up and doubles the total correction power.
        '''
        # erase_pos is a list of positions where you know (or greatly suspect at least) there is an erasure (ie, wrong character but you know it's at this position). Just input the list of all positions you know there are errors, and this method will automatically split the erasures positions to attach to the corresponding data chunk.

        cdef ca.array data_array = <ca.array>pa.array('Q', data)
        cdef uint64_t[::1] data_array_view = data_array

        cdef bool has_erase_pos_array = False
        cdef ca.array erase_pos_array
        if erase_pos is not None:
            has_erase_pos_array = True
            erase_pos_array = <ca.array>pa.array('Q', erase_pos)

        # Precache class attributes into local variables
        cdef uint64_t nsym_c
        if nsym is None:
            nsym_c = self.nsym
        else:
            nsym_c = nsym

        cdef uint64_t nsize = self.nsize
        cdef uint64_t fcr = self.fcr
        cdef uint64_t generator = self.generator

        # Calculate chunk size and total number of chunks for looping
        cdef uint64_t chunk_size = nsize
        cdef uint64_t total_chunks = (<uint64_t>len(data_array) - 1) // chunk_size + 1  # need to convert to floats first to get an accurate floating division, or else we assume implicit conversion and it will cause an error on Python 2
        cdef uint64_t nmes = nsize - nsym_c

        # Initialize output array
        cdef ca.array dec = ca.clone(uint64_array_template, <Py_ssize_t>(nmes * total_chunks), False)  # pre-allocate array and we will overwrite data in it, much faster than extending
        cdef uint64_t[::1] dec_view = dec
        cdef ca.array dec_full = ca.clone(uint64_array_template, <Py_ssize_t>(nsize * total_chunks), False)
        cdef uint64_t[::1] dec_full_view = dec_full
        cdef ca.array errata_pos_all = ca.clone(uint64_array_template, <Py_ssize_t>(nsize * total_chunks), False)
        cdef uint64_t errata_pos_all_view_size = 0

        # Chunking loop
        cdef uint64_t i , j
        cdef ca.array e_pos, tmp, rmes, recc, errata_pos
        cdef uint64_t e_pos_size, tmp_size
        cdef uint64_t[::1] rmes_view, recc_view
        for i in range(total_chunks):  # Split the long message in a chunk
            # Extract the erasures for this chunk
            if has_erase_pos_array:
                # First extract the erasures for this chunk (all erasures below the maximum chunk length)
                e_pos = ca.clone(uint64_array_template, len(erase_pos_array), False)
                e_pos_size = 0
                for j in range(<uint64_t>len(erase_pos_array)):
                    if erase_pos_array.data.as_ulonglongs[j] < nsize:
                        e_pos.data.as_ulonglongs[e_pos_size] = erase_pos_array.data.as_ulonglongs[j]
                        e_pos_size += 1
                ca.resize(e_pos, e_pos_size)
                # Then remove the extract erasures from the big list and also decrement all subsequent positions values by nsize (the current chunk's size) so as to prepare the correct alignment for the next iteration
                tmp = ca.clone(uint64_array_template, len(erase_pos_array), False)
                tmp_size = 0
                for j in range(<uint64_t>len(erase_pos_array)):
                    if erase_pos_array.data.as_ulonglongs[j] >= nsize:
                        tmp.data.as_ulonglongs[tmp_size] = erase_pos_array.data.as_ulonglongs[j] - nsize
                        tmp_size += 1
                erase_pos_array = tmp
            else:
                e_pos = ca.clone(uint64_array_template, 0, False)
            # Decode/repair this chunk!
            rmes, recc, errata_pos = rs_correct_msg(self.gf, data_array_view[chunk_size * i:chunk_size * (i + 1)], nsym_c, fcr, generator, e_pos, only_erasures)
            rmes_view, recc_view = rmes, recc
            dec_view[nmes * i:nmes * i + rmes_view.shape[0]] = rmes_view
            if i == total_chunks - 1:
                ca.resize(dec, <Py_ssize_t>(nmes * i + <uint64_t>rmes_view.shape[0]))
            dec_full_view[nsize * i:nsize * i + rmes_view.shape[0]] = rmes_view
            if i == total_chunks - 1:
                ca.resize(dec_full, <Py_ssize_t>(nsize * i + <uint64_t>rmes_view.shape[0] + <uint64_t>recc_view.shape[0]))
                dec_full_view = dec_full
            dec_full_view[nsize * i + rmes_view.shape[0]:nsize * i + rmes_view.shape[0] + recc_view.shape[0]] = recc_view  # append corrected ecc just after corrected message. The two lines are equivalent to rmes + recc but here we don't need to concatenate both arrays first (and create a third one for nothing) before storing in the output array
            for j in range(<uint64_t>len(errata_pos)):
                errata_pos_all.data.as_ulonglongs[errata_pos_all_view_size] = errata_pos.data.as_ulonglongs[j] + i * nsize
                errata_pos_all_view_size += 1
        ca.resize(errata_pos_all, errata_pos_all_view_size)
        return dec, dec_full, errata_pos_all

    def check(self, data: Iterable[int], nsym: Optional[int] = None) -> Annotated[pa.array[int], 'B']:
        '''Check if a message+ecc stream is not corrupted (or fully repaired). Note: may return a wrong result if number of errors > nsym.'''

        cdef uint64_t nsym_c
        if nsym is None:
            nsym_c = self.nsym
        else:
            nsym_c = nsym

        cdef ca.array data_array = <ca.array>pa.array('Q', data)
        cdef uint64_t[::1] data_array_view = data_array

        # Precache class attributes into local variables
        cdef uint64_t nsize = self.nsize
        cdef uint64_t fcr = self.fcr
        cdef uint64_t generator = self.generator

        # Calculate chunksize
        cdef uint64_t chunk_size = nsize
        cdef uint64_t total_chunks = (<uint64_t>len(data_array) - 1) // chunk_size + 1  # need to convert to floats first to get an accurate floating division, or else we assume implicit conversion and it will cause an error on Python 2

        # Pre-allocate output array
        cdef ca.array check = ca.clone(uint8_array_template, total_chunks, False)  # TODO: could use a bitarray but this creates an external dependency and it's not compatible with cython https://pypi.org/project/bitarray/ and https://www.noveltech.dev/booleans-python-numpy/

        # Chunking loop
        cdef uint64_t i
        for i in range(total_chunks):  # Split the long message in a chunk
            check.data.as_uchars[i] = rs_check(self.gf, data_array_view[chunk_size * i:chunk_size * (i + 1)], nsym_c, fcr, generator)
        return check

    cpdef tuple maxerrata(self, uint64_t nsym = <uint64_t>(-1), uint64_t errors = <uint64_t>(-1), uint64_t erasures = <uint64_t>(-1), bool verbose = False):
        '''Return the Singleton Bound for the current codec, which is the max number of errata (errors and erasures) that the codec can decode/correct.
        Beyond the Singleton Bound (too many errors/erasures), the algorithm will try to raise an exception, but it may also not detect any problem with the message and return 0 errors.
        Hence why you should use checksums if your goal is to detect errors (as opposed to correcting them), as checksums have no bounds on the number of errors, the only limitation being the probability of collisions.
        By default, return a tuple wth the maximum number of errors (2nd output) OR erasures (2nd output) that can be corrected.
        If errors or erasures (not both) is specified as argument, computes the remaining **simultaneous** correction capacity (eg, if errors specified, compute the number of erasures that can be simultaneously corrected).
        Set verbose to True to get print a report.'''
        # Fetch nsym from class attributes if not overriden by a function call
        cdef uint64_t nsym_c
        if nsym != <uint64_t>(-1):
            nsym_c = self.nsym
        else:
            nsym_c = nsym
        # Compute the maximum number of errors OR erasures
        cdef uint64_t maxerrors = nsym_c // 2  # always floor the number, we can't correct half a symbol, it's all or nothing
        cdef uint64_t maxerasures = nsym_c
        # Compute the maximum of simultaneous errors AND erasures
        if erasures != <uint64_t>(-1):
            # We know the erasures count, we want to know how many errors we can correct simultaneously
            if erasures > maxerasures:
                raise ReedSolomonError('Specified number of errors or erasures exceeding the Singleton Bound!')
            maxerrors = (nsym_c - erasures) // 2
            if verbose:
                print(f'This codec can correct up to {maxerrors} errors and {erasures} erasures simultaneously')
            # Return a tuple with the maximum number of simultaneously corrected errors and erasures
            return maxerrors, erasures
        if errors != <uint64_t>(-1):
            # We know the errors count, we want to know how many erasures we can correct simultaneously
            if errors > maxerrors:
                raise ReedSolomonError('Specified number of errors or erasures exceeding the Singleton Bound!')
            maxerasures = nsym_c - (errors * 2)
            if verbose:
                print(f'This codec can correct up to {errors} errors and {maxerasures} erasures simultaneously')
            # Return a tuple with the maximum number of simultaneously corrected errors and erasures
            return errors, maxerasures
        # Return a tuple with the maximum number of errors and erasures (independently corrected)
        if verbose:
            print(f'This codec can correct up to {maxerrors} errors and {maxerasures} erasures independently')
        return maxerrors, maxerasures
