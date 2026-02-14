cimport cpython.array as ca
from libc.stdint cimport *

################### INIT and stuff ###################

cdef ca.array uint64_array_template
cdef ca.array uint8_array_template


################### GALOIS FIELD ELEMENTS MATHS ###################

cpdef ca.array rwh_primes1(uint64_t n)

# cpdef uint64_t cl_mult(uint64_t x, uint64_t y)
# cpdef uint64_t bit_length(uint64_t n)
# cpdef uint64_t cl_div(uint64_t dividend, uint64_t divisor)
# cpdef uint64_t gf_mult_noLUT_slow(uint64_t x, uint64_t y, uint64_t prim = *)

cpdef uint64_t gf_mult_noLUT(uint64_t x, uint64_t y, uint64_t prim = *, uint64_t field_charac_full = *, bool carryless = *)

cpdef ca.array find_prime_polys(uint64_t generator = *, uint64_t c_exp = *, bool fast_primes = *, bool single = *)

cdef class GaloisField:
    cdef uint64_t field_charac
    cdef ca.array gf_exp, gf_log
    cpdef uint64_t gf_inverse(self, uint64_t x)
    cpdef uint64_t gf_mul(self, uint64_t x, uint64_t y)
    cpdef uint64_t gf_div(self, uint64_t x, uint64_t y)
    cpdef uint64_t gf_pow(self, uint64_t x, int64_t power)

cpdef uint64_t gf_add(uint64_t x, uint64_t y)
cpdef uint64_t gf_sub(uint64_t x, uint64_t y)
cpdef uint64_t gf_neg(uint64_t x)


################### GALOIS FIELD POLYNOMIALS MATHS ###################

cpdef ca.array gf_poly_scale(GaloisField gf, const uint64_t[::1] p, uint64_t x)
cpdef ca.array gf_poly_add(const uint64_t[::1] p, const uint64_t[::1] q)
cpdef ca.array gf_poly_mul(GaloisField gf, const uint64_t[::1] p, const uint64_t[::1] q)
# cpdef ca.array gf_poly_mul_simple(GaloisField gf, const uint64_t[::1] p, const uint64_t[::1] q)
cpdef ca.array gf_poly_neg(const uint64_t[::1] poly)
cpdef tuple gf_poly_div(GaloisField gf, const uint64_t[::1] dividend, const uint64_t[::1] divisor)
cpdef ca.array gf_poly_square(GaloisField gf, const uint64_t[::1] poly)
cpdef uint64_t gf_poly_eval(GaloisField gf, const uint64_t[::1] poly, uint64_t x)


################### REED-SOLOMON ENCODING ###################

cpdef ca.array rs_generator_poly(GaloisField gf, uint64_t nsym, uint64_t fcr, uint64_t generator = *)
cpdef list rs_generator_poly_all(GaloisField gf, uint64_t max_nsym, uint64_t fcr, uint64_t generator = *)
# cpdef ca.array rs_simple_encode_msg(GaloisField gf, const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *, const uint64_t[::1] gen = *)
cpdef ca.array rs_encode_msg(GaloisField gf, const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *, const uint64_t[::1] gen = *)


################### REED-SOLOMON DECODING ###################

cpdef ca.array inverted(const uint64_t[::1] msg)
cpdef ca.array rs_calc_syndromes(GaloisField gf, const uint64_t[::1] msg, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *)
cpdef ca.array rs_correct_errata(GaloisField gf, const uint64_t[::1] msg_in, const uint64_t[::1] synd, const uint64_t[::1] err_pos, uint64_t fcr = *, uint64_t generator = *)
cpdef ca.array rs_find_error_locator(GaloisField gf, const uint64_t[::1] synd, uint64_t nsym, const uint64_t[::1] erase_loc = *, uint64_t erase_count = *)
cpdef ca.array rs_find_errata_locator(GaloisField gf, const uint64_t[::1] e_pos, uint64_t generator = *)
cpdef ca.array rs_find_error_evaluator(GaloisField gf, const uint64_t[::1] synd, const uint64_t[::1] err_loc, uint64_t nsym)
cpdef ca.array rs_find_errors(GaloisField gf, const uint64_t[::1] err_loc, uint64_t nmess, uint64_t generator = *)
cpdef ca.array rs_forney_syndromes(GaloisField gf, const uint64_t[::1] synd, const uint64_t[::1] pos, uint64_t nmess, uint64_t generator = *)
cpdef tuple rs_correct_msg(GaloisField gf, const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *, const uint64_t[::1] erase_pos = *, bool only_erasures = *)
# cpdef tuple rs_correct_msg_nofsynd(const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *, const uint64_t[::1] erase_pos = *, bool only_erasures = *)
cpdef bool rs_check(GaloisField gf, const uint64_t[::1] msg, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *)


#===================================================================================================
# API
#===================================================================================================
cdef class RSCodec:
    cdef uint64_t nsym, nsize, fcr, prim, generator, c_exp
    cdef GaloisField gf
    cdef bool single_gen
    cdef ca.array gen
    cdef list gens

    cpdef tuple maxerrata(self, uint64_t nsym = *, uint64_t errors = *, uint64_t erasures = *, bool verbose = *)
