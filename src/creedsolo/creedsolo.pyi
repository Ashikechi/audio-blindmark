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
    >> rmes, rmesecc = rsc.decode(mesecc)
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

import array
from typing import Annotated, Generator, Iterable, Optional, Sequence

################### INIT and stuff ###################

class ReedSolomonError(Exception):
    pass

################### GALOIS FIELD ELEMENTS MATHS ###################
# General note: Galois Field maths essentially are all the standard math operations everybody learn in primary school,
# but with only integer AND they can wraparound (ie, we use a modulo), so that in practice this means that
# Galois Field math operations are bounded in a very specific range of values. This changes a lot how the maths are done,
# but not that much, so you can still wrap your head around if you are willing to spend some time.

def rwh_primes1(n: int) -> Annotated[array.array, 'Q']:
    ''' Returns a list of primes < n '''
    ...
# def gf_mult_noLUT_slow(x: int, y: int, prim: Optional[int] = None) -> int:  # pylint: disable=C0103
#     '''Multiplication in Galois Fields on-the-fly without using a precomputed look-up table (and thus it's slower) by using the standard carry-less multiplication + modular reduction using an irreducible prime polynomial.'''
#     ...
def gf_mult_noLUT(x: int, y: int, prim: int = 0, field_charac_full: int = 256, carryless: bool = True) -> int:  # pylint: disable=C0103
    '''Galois Field integer multiplication on-the-fly without using a look-up table, using Russian Peasant Multiplication algorithm (faster than the standard multiplication + modular reduction). This is still slower than using a look-up table, but is the fastest alternative, and is often used in embedded circuits where storage space is limited (ie, no space for a look-up table).
    If prim is 0 and carryless=False, then the function produces the result for a standard integers multiplication (no carry-less arithmetics nor modular reduction).'''
    ...
def find_prime_polys(generator: int = 2, c_exp: int = 8, fast_primes: bool = False, single: bool = False) -> Annotated[array.array, 'Q']:
    '''Compute the list of prime polynomials for the given generator and galois field characteristic exponent.'''
    ...

class GaloisField:
    def __init__(self, prim: int = 0x11d, generator: int = 2, c_exp: int = 8) -> None:
        '''Precompute the logarithm and anti-log tables for faster computation later, using the provided primitive polynomial.
        These tables are used for multiplication/division since addition/substraction are simple XOR operations inside GF of characteristic 2.
        The basic idea is quite simple: since b ** (log_b(x), log_b(y)) == x * y given any number b (the base or generator of the logarithm), then we can use any number b to precompute logarithm and anti-log (exponentiation) tables to use for multiplying two numbers x and y.
        That's why when we use a different base/generator number, the log and anti-log tables are drastically different, but the resulting computations are the same given any such tables.
        For more infos, see https://en.wikipedia.org/wiki/Finite_field_arithmetic#Implementation_tricks
        '''
        ...
    def gf_inverse(self, x: int) -> int:
        '''Inverse of a galois field integer'''
        ...
    def gf_mul(self, x: int, y: int) -> int:
        '''Multiply two galois field integers'''
        ...
    def gf_div(self, x: int, y: int) -> int:
        '''Divide x by y galois field integers'''
        ...
    def gf_pow(self, x: int, power: int) -> int:
        '''Power of x galois field integer'''
        ...

def gf_add(x: int, y: int) -> int:
    '''Add two galois field integers'''
    ...
def gf_sub(x: int, y: int) -> int:
    '''Subtract two galois field integers'''
    ...
def gf_neg(x: int) -> int:
    '''Negate one galois field integer (does nothing)'''
    ...

################### GALOIS FIELD POLYNOMIALS MATHS ###################

def gf_poly_scale(gf: GaloisField, p: Annotated[array.array, 'Q'], x: int) -> Annotated[array.array, 'Q']:
    '''Scale a galois field polynomial with a factor x (an integer)'''
    ...
def gf_poly_add(p: Annotated[array.array, 'Q'], q: Annotated[array.array, 'Q']) -> Annotated[array.array, 'Q']:
    '''Add two galois field polynomials'''
    ...
def gf_poly_mul(gf: GaloisField, p: Annotated[array.array, 'Q'], q: Annotated[array.array, 'Q']) -> Annotated[array.array, 'Q']:
    '''Multiply two polynomials, inside Galois Field (but the procedure is generic). Optimized function by precomputation of log.'''
    ...
# def gf_poly_mul_simple(gf: GaloisField, p: Annotated[array.array, 'Q'], q: Annotated[array.array, 'Q']): # simple equivalent way of multiplying two polynomials without precomputation, but thus it's slower
#     '''Multiply two polynomials, inside Galois Field'''
#     ...
def gf_poly_neg(poly: Annotated[array.array, 'Q']) -> Annotated[array.array, 'Q']:
    '''Returns the polynomial with all coefficients negated. In GF(2^p), negation does not change the coefficient, so we return the polynomial as-is.'''
    ...
def gf_poly_div(gf: GaloisField, dividend: Annotated[array.array, 'Q'], divisor: Annotated[array.array, 'Q']) -> tuple[Annotated[array.array, 'Q'], Annotated[array.array, 'Q']]:
    '''Fast polynomial division by using Extended Synthetic Division and optimized for GF(2^p) computations (doesn't work with standard polynomials outside of this galois field).'''
    ...
def gf_poly_square(gf: GaloisField, poly: Annotated[array.array, 'Q']) -> Annotated[array.array, 'Q']:  # pragma: no cover
    '''Linear time implementation of polynomial squaring. For details, see paper: "A fast software implementation for arithmetic operations in GF (2n)". De Win, E., Bosselaers, A., Vandenberghe, S., De Gersem, P., & Vandewalle, J. (1996, January). In Advances in Cryptology - Asiacrypt'96 (pp. 65-76). Springer Berlin Heidelberg.'''
    ...
def gf_poly_eval(gf: GaloisField, poly: Annotated[array.array, 'Q'], x: int) -> int:
    '''Evaluates a polynomial in GF(2^p) given the value for x. This is based on Horner's scheme for maximum efficiency.'''
    ...

################### REED-SOLOMON ENCODING ###################

def rs_generator_poly(gf: GaloisField, nsym: int, fcr: int = 0, generator: int = 2) -> Annotated[array.array, 'Q']:
    '''Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)'''
    ...
def rs_generator_poly_all(gf: GaloisField, max_nsym: int, fcr: int = 0, generator: int = 2) -> list[Annotated[array.array, 'Q']]:
    '''Generate all irreducible generator polynomials up to max_nsym (usually you can use n, the length of the message+ecc). Very useful to reduce processing time if you want to encode using variable schemes and nsym rates.'''
    ...
# def rs_simple_encode_msg(gf: GaloisField, msg_in: Annotated[array.array, 'Q'], nsym: int, fcr: int = 0, generator: int = 2, gen: Optional[Annotated[array.array, 'Q']] = None) -> Annotated[array.array, 'Q']:
#     '''Simple Reed-Solomon encoding (mainly an example for you to understand how it works, because it's slower than the inlined function below)'''
#     ...
def rs_encode_msg(gf: GaloisField, msg_in: Annotated[array.array, 'Q'], nsym: int, fcr: int = 0, generator: int = 2, gen: Optional[Annotated[array.array, 'Q']] = None) -> Annotated[array.array, 'Q']:
    '''Reed-Solomon main encoding function, using polynomial division (Extended Synthetic Division, the fastest algorithm available to my knowledge), better explained at http://research.swtch.com/field'''
    ...

################### REED-SOLOMON DECODING ###################

def inverted(msg: Annotated[array.array, 'Q']) -> Annotated[array.array, 'Q']:
    '''Implements msg[::-1] explicitly to make the library compatible with MicroPython which does not support stepped slices.'''
    ...
def rs_calc_syndromes(gf: GaloisField, msg: Annotated[array.array, 'Q'], nsym: int, fcr: int = 0, generator: int = 2) -> Annotated[array.array, 'Q']:
    '''Given the received codeword msg and the number of error correcting symbols (nsym), computes the syndromes polynomial.
    Mathematically, it's essentially equivalent to a Fourrier Transform (Chien search being the inverse).
    '''
    ...
def rs_correct_errata(gf: GaloisField, msg_in: Annotated[array.array, 'Q'], synd: Annotated[array.array, 'Q'], err_pos: Annotated[array.array, 'Q'], fcr: int = 0, generator: int = 2) -> Annotated[array.array, 'Q']: # err_pos is a list of the positions of the errors/erasures/errata
    '''Forney algorithm, computes the values (error magnitude) to correct the input message.'''
    ...
def rs_find_error_locator(gf: GaloisField, synd: Annotated[array.array, 'Q'], nsym: int, erase_loc: Optional[Annotated[array.array, 'Q']] = None, erase_count: int = 0) -> Annotated[array.array, 'Q']:
    '''Find error/errata locator and evaluator polynomials with Berlekamp-Massey algorithm'''
    ...
def rs_find_errata_locator(gf: GaloisField, e_pos: Annotated[array.array, 'Q'], generator: int = 2) -> Annotated[array.array, 'Q']:
    '''Compute the erasures/errors/errata locator polynomial from the erasures/errors/errata positions (the positions must be relative to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to "h_ll_ worldxxxxxxxxx" with xxxxxxxxx being the ecc of length n-k=9, here the string positions are [1, 4], but the coefficients are reversed since the ecc characters are placed as the first coefficients of the polynomial, thus the coefficients of the erased characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified as an argument.'''
    ...
def rs_find_error_evaluator(gf: GaloisField, synd: Annotated[array.array, 'Q'], err_loc: Annotated[array.array, 'Q'], nsym: int) -> Annotated[array.array, 'Q']:
    '''Compute the error (or erasures if you supply sigma=erasures locator polynomial, or errata) evaluator polynomial Omega from the syndrome and the error/erasures/errata locator Sigma. Omega is already computed at the same time as Sigma inside the Berlekamp-Massey implemented above, but in case you modify Sigma, you can recompute Omega afterwards using this method, or just ensure that Omega computed by BM is correct given Sigma.'''
    ...
def rs_find_errors(gf: GaloisField, err_loc: Annotated[array.array, 'Q'], nmess: int, generator: int = 2) -> Annotated[array.array, 'Q']:
    '''Find the roots (ie, where evaluation = zero) of error polynomial by smart bruteforce trial. This is a faster form of chien search, processing only useful coefficients (the ones in the messages) instead of the whole 2^8 range. Besides the speed boost, this also allows to fix a number of issue: correctly decoding when the last ecc byte is corrupted, and accepting messages of length n > 2^8.'''
    ...
def rs_forney_syndromes(gf: GaloisField, synd: Annotated[array.array, 'Q'], pos: Annotated[array.array, 'Q'], nmess: int, generator: int = 2) -> Annotated[array.array, 'Q']:
    ...
def rs_correct_msg(gf: GaloisField, msg_in: Annotated[array.array, 'Q'], nsym: int, fcr: int = 0, generator: int = 2, erase_pos: Optional[Annotated[array.array, 'Q']] = None, only_erasures: bool = False) -> tuple[Annotated[array.array, 'Q'], Annotated[array.array, 'Q'], Annotated[array.array, 'Q']]:
    '''Reed-Solomon main decoding function'''
    ...
# def rs_correct_msg_nofsynd(msg_in: Annotated[array.array, 'Q'], nsym: int, fcr: int = 0, generator: int = 2, erase_pos: Optional[Annotated[array.array, 'Q']] = None, only_erasures: bool = False) -> tuple[Annotated[array.array, 'Q'], Annotated[array.array, 'Q'], Annotated[array.array, 'Q']]:
#     '''Reed-Solomon main decoding function, without using the modified Forney syndromes'''
#     ...
def rs_check(gf: GaloisField, msg: Annotated[array.array, 'Q'], nsym: int, fcr: int = 0, generator: int = 2) -> bool:
    '''Returns true if the message + ecc has no error of false otherwise (may not always catch a wrong decoding or a wrong message, particularly if there are too many errors -- above the Singleton bound --, but it usually does)'''
    ...

######################## end of REED SOLOMON DECODING ###############

#===================================================================================================
# API
#===================================================================================================
class RSCodec:
    '''
    A Reed Solomon encoder/decoder. After initializing the object, use ``encode`` to encode a
    (byte)string to include the RS correction code, and pass such an encoded (byte)string to
    ``decode`` to extract the original message (if the number of errors allows for correct decoding).
    The ``nsym`` argument is the length of the correction code, and it determines the number of
    error bytes (if I understand this correctly, half of ``nsym`` is correctable)

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
        ...
    def chunk(self, data: Sequence[int], chunk_size: int) -> Generator[Sequence[int], None, None]:
        '''Split a long message into chunks
        DEPRECATED: inlined alternate form so that we can preallocate arrays and hence get faster results with JIT compilers such as PyPy.'''
        ...
    def encode(self, data: Iterable[int], nsym: Optional[int] = None) -> Annotated[array.array, 'Q']:
        '''Encode a message (ie, add the ecc symbols) using Reed-Solomon, whatever the length of the message because we use chunking
        Optionally, can set nsym to encode with a different number of error correction symbols, but RSCodec must be initialized with single_gen=False first.
        slice_assign=True allows to speed up the loop quite significantly in JIT compilers such as PyPy by preallocating the output bytearray and slice assigning into it, instead of constantly extending an empty bytearray, but this only works in Python 3, not Python 2, hence is disabled by default for retrocompatibility.
        '''
        ...
    def decode(self, data: Iterable[int], nsym: Optional[int] = None, erase_pos: Optional[Iterable[int]] = None, only_erasures: bool = False) -> tuple[Annotated[array.array, 'Q'], Annotated[array.array, 'Q'], Annotated[array.array, 'Q']]:
        '''Repair a message, whatever its size is, by using chunking. May return a wrong result if number of errors > nsym because then too many errors to be corrected.
        Note that it returns a couple of vars: the repaired messages, and the repaired messages+ecc (useful for checking).
        Usage: rmes, rmesecc = RSCodec.decode(data).
        Optionally: can specify nsym to decode messages of different parameters, erase_pos with a list of erasures positions to double the number of erasures that can be corrected compared to unlocalized errors, only_erasures boolean to specify if we should only look for erasures, which speeds up and doubles the total correction power.
        '''
        ...
    def check(self, data: Iterable[int], nsym: Optional[int] = None) -> Annotated[array.array, 'B']:
        '''Check if a message+ecc stream is not corrupted (or fully repaired). Note: may return a wrong result if number of errors > nsym.'''
        ...
    def maxerrata(self, nsym: Optional[int] = None, errors: Optional[int] = None, erasures: Optional[int] = None, verbose: bool = False) -> tuple[int, int]:
        '''Return the Singleton Bound for the current codec, which is the max number of errata (errors and erasures) that the codec can decode/correct.
        Beyond the Singleton Bound (too many errors/erasures), the algorithm will try to raise an exception, but it may also not detect any problem with the message and return 0 errors.
        Hence why you should use checksums if your goal is to detect errors (as opposed to correcting them), as checksums have no bounds on the number of errors, the only limitation being the probability of collisions.
        By default, return a tuple wth the maximum number of errors (2nd output) OR erasures (2nd output) that can be corrected.
        If errors or erasures (not both) is specified as argument, computes the remaining **simultaneous** correction capacity (eg, if errors specified, compute the number of erasures that can be simultaneously corrected).
        Set verbose to True to get print a report.'''
        ...
