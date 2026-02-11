from libc.stdint cimport * # pyright: ignore[reportWildcardImportFromLibrary]
cimport cpython.array as ca


################### INIT and stuff ###################

cdef ca.array uint64_array_template
cdef ca.array uint8_array_template

cdef ca.array gf_exp
cdef ca.array gf_log
cdef uint64_t field_charac


################### GALOIS FIELD ELEMENTS MATHS ###################

cpdef ca.array rwh_primes1(uint64_t n)
cpdef ca.array find_prime_polys(uint64_t generator = *, uint64_t c_exp = *, bool fast_primes = *, bool single = *)
cpdef tuple init_tables(uint64_t prim = *, uint64_t generator = *, uint64_t c_exp = *)

cpdef uint64_t gf_add(uint64_t x, uint64_t y)
cpdef uint64_t gf_sub(uint64_t x, uint64_t y)
cpdef uint64_t gf_neg(uint64_t x)
cpdef uint64_t gf_inverse(uint64_t x)
cpdef uint64_t gf_mul(uint64_t x, uint64_t y)
cpdef uint64_t gf_div(uint64_t x, uint64_t y)
cpdef uint64_t gf_pow(uint64_t x, int64_t power)

# cpdef uint64_t cl_mult(uint64_t x, uint64_t y)
# cpdef uint64_t bit_length(uint64_t n)
# cpdef uint64_t cl_div(uint64_t dividend, uint64_t divisor)
# cpdef uint64_t gf_mult_noLUT_slow(uint64_t x, uint64_t y, uint64_t prim = *)

cpdef uint64_t gf_mult_noLUT(uint64_t x, uint64_t y, uint64_t prim = *, uint64_t field_charac_full = *, bool carryless = *)


################### GALOIS FIELD POLYNOMIALS MATHS ###################

cpdef ca.array gf_poly_scale(const uint64_t[::1] p, uint64_t x)
cpdef ca.array gf_poly_add(const uint64_t[::1] p, const uint64_t[::1] q)
cpdef ca.array gf_poly_mul(const uint64_t[::1] p, const uint64_t[::1] q)
# cpdef ca.array gf_poly_mul_simple(const uint64_t[::1] p, const uint64_t[::1] q)
cpdef ca.array gf_poly_neg(const uint64_t[::1] poly)
cpdef tuple gf_poly_div(const uint64_t[::1] dividend, const uint64_t[::1] divisor)
cpdef uint64_t gf_poly_eval(const uint64_t[::1] poly, uint64_t x)


################### REED-SOLOMON ENCODING ###################

cpdef ca.array rs_generator_poly(uint64_t nsym, uint64_t fcr, uint64_t generator = *)
cpdef list rs_generator_poly_all(uint64_t max_nsym, uint64_t fcr, uint64_t generator = *)
# cpdef ca.array rs_simple_encode_msg(const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *, const uint64_t[::1] gen = *)
cpdef ca.array rs_encode_msg(const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *, const uint64_t[::1] gen = *)


################### REED-SOLOMON DECODING ###################

cpdef ca.array inverted(const uint64_t[::1] msg)
cpdef ca.array rs_calc_syndromes(const uint64_t[::1] msg, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *)
cpdef ca.array rs_correct_errata(const uint64_t[::1] msg_in, const uint64_t[::1] synd, const uint64_t[::1] err_pos, uint64_t fcr = *, uint64_t generator = *)
cpdef ca.array rs_find_error_locator(const uint64_t[::1] synd, uint64_t nsym, const uint64_t[::1] erase_loc = *, uint64_t erase_count = *)
cpdef ca.array rs_find_errata_locator(const uint64_t[::1] e_pos, uint64_t generator = *)
cpdef ca.array rs_find_error_evaluator(const uint64_t[::1] synd, const uint64_t[::1] err_loc, uint64_t nsym)
cpdef ca.array rs_find_errors(const uint64_t[::1] err_loc, uint64_t nmess, uint64_t generator = *)
cpdef ca.array rs_forney_syndromes(const uint64_t[::1] synd, const uint64_t[::1] pos, uint64_t nmess, uint64_t generator = *)
cpdef tuple rs_correct_msg(const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *, const uint64_t[::1] erase_pos = *, bool only_erasures = *)
# cpdef tuple rs_correct_msg_nofsynd(const uint64_t[::1] msg_in, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *, const uint64_t[::1] erase_pos = *, bool only_erasures = *)
cpdef bool rs_check(const uint64_t[::1] msg, uint64_t nsym, uint64_t fcr = *, uint64_t generator = *)


#===================================================================================================
# API
#===================================================================================================
cdef class RSCodec:
    cdef uint64_t nsym, nsize, fcr, prim, generator, c_exp, field_charac
    cdef ca.array gf_log, gf_exp
    cdef bool single_gen
    cdef ca.array gen
    cdef list gens

    cdef object __globals_temp

    cpdef tuple maxerrata(self, uint64_t nsym = *, uint64_t errors = *, uint64_t erasures = *, bool verbose = *) # pyright: ignore[reportGeneralTypeIssues]
