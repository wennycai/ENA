/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (C) 2007-2019 RELIC Authors
 *
 * This file is part of RELIC. RELIC is legal property of its developers,
 * whose names are not listed here. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * RELIC is free software; you can redistribute it and/or modify it under the
 * terms of the version 2.1 (or later) of the GNU Lesser General Public License
 * as published by the Free Software Foundation; or version 2.0 of the Apache
 * License as published by the Apache Software Foundation. See the LICENSE files
 * for more details.
 *
 * RELIC is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the LICENSE files for more details.
 *
 * You should have received a copy of the GNU Lesser General Public or the
 * Apache License along with RELIC. If not, see <https://www.gnu.org/licenses/>
 * or <https://www.apache.org/licenses/>.
 */

/**
 * @file
 *
 * Tests for extensions defined over prime fields.
 *
 * @ingroup test
 */

#include "relic.h"
#include "relic_test.h"

static int memory2(void) {
	err_t e;
	int code = RLC_ERR;
	fp2_t a;

	fp2_null(a);

	TRY {
		TEST_BEGIN("memory can be allocated") {
			fp2_new(a);
			fp2_free(a);
		} TEST_END;
	} CATCH(e) {
		switch (e) {
			case ERR_NO_MEMORY:
				util_print("FATAL ERROR!\n");
				ERROR(end);
				break;
		}
	}
	(void)a;
	code = RLC_OK;
  end:
	return code;
}

static int util2(void) {
	int code = RLC_ERR;
	uint8_t bin[2 * RLC_FP_BYTES];
	fp2_t a, b, c;
	dig_t d;

	fp2_null(a);
	fp2_null(b);
	fp2_null(c);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);

		TEST_BEGIN("comparison is consistent") {
			fp2_rand(a);
			fp2_rand(b);
			if (fp2_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp2_cmp(b, a) == RLC_NE, end);
			}
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_rand(c);
			if (fp2_cmp(a, c) != RLC_EQ) {
				fp2_copy(c, a);
				TEST_ASSERT(fp2_cmp(c, a) == RLC_EQ, end);
			}
			if (fp2_cmp(b, c) != RLC_EQ) {
				fp2_copy(c, b);
				TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
			}
		}
		TEST_END;

		TEST_BEGIN("negation is consistent") {
			fp2_rand(a);
			fp2_neg(b, a);
			if (fp2_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp2_cmp(b, a) == RLC_NE, end);
			}
			fp2_neg(b, b);
			TEST_ASSERT(fp2_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and comparison are consistent") {
			do {
				fp2_rand(a);
			} while (fp2_is_zero(a));
			fp2_zero(c);
			TEST_ASSERT(fp2_cmp(a, c) == RLC_NE, end);
			TEST_ASSERT(fp2_cmp(c, a) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to random and comparison are consistent") {
			do {
				fp2_rand(a);
			} while (fp2_is_zero(a));
			fp2_rand(c);
			TEST_ASSERT(fp2_cmp(a, c) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and zero test are consistent") {
			fp2_zero(a);
			TEST_ASSERT(fp2_is_zero(a), end);
		}
		TEST_END;

		TEST_BEGIN("assignment to a constant and comparison are consistent") {
			rand_bytes((uint8_t *)&d, (RLC_DIG / 8));
			fp2_set_dig(a, d);
			TEST_ASSERT(fp2_cmp_dig(a, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("reading and writing a finite field element are consistent") {
			fp2_rand(a);
			fp2_write_bin(bin, sizeof(bin), a, 0);
			fp2_read_bin(b, bin, sizeof(bin));
			TEST_ASSERT(fp2_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("getting the size of a prime field element is correct") {
			fp2_rand(a);
			TEST_ASSERT(fp2_size_bin(a, 0) == 2 * RLC_FP_BYTES, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	return code;
}

static int addition2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c, d, e;

	fp2_null(a);
	fp2_null(b);
	fp2_null(c);
	fp2_null(d);
	fp2_null(e);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);
		fp2_new(d);
		fp2_new(e);

		TEST_BEGIN("addition is commutative") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_add(d, a, b);
			fp2_add(e, b, a);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition is associative") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_rand(c);
			fp2_add(d, a, b);
			fp2_add(d, d, c);
			fp2_add(e, b, c);
			fp2_add(e, a, e);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has identity") {
			fp2_rand(a);
			fp2_zero(d);
			fp2_add(e, a, d);
			TEST_ASSERT(fp2_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has inverse") {
			fp2_rand(a);
			fp2_neg(d, a);
			fp2_add(e, a, d);
			TEST_ASSERT(fp2_is_zero(e), end);
		} TEST_END;

#if PP_QDR == BASIC || !defined(STRIP)
		TEST_BEGIN("basic addition is correct") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_add(d, a, b);
			fp2_add_basic(e, a, b);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_QDR == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated addition is correct") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_add(d, a, b);
			fp2_add_integ(e, a, b);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	fp2_free(d);
	fp2_free(e);
	return code;
}

static int subtraction2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c, d;

	fp2_null(a);
	fp2_null(b);
	fp2_null(c);
	fp2_null(d);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);
		fp2_new(d);

		TEST_BEGIN("subtraction is anti-commutative") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_sub(c, a, b);
			fp2_sub(d, b, a);
			fp2_neg(d, d);
			TEST_ASSERT(fp2_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has identity") {
			fp2_rand(a);
			fp2_zero(c);
			fp2_sub(d, a, c);
			TEST_ASSERT(fp2_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has inverse") {
			fp2_rand(a);
			fp2_sub(c, a, a);
			TEST_ASSERT(fp2_is_zero(c), end);
		}
		TEST_END;

#if PP_QDR == BASIC || !defined(STRIP)
		TEST_BEGIN("basic subtraction is correct") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_sub(c, a, b);
			fp2_sub_basic(d, a, b);
			TEST_ASSERT(fp2_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_QDR == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated subtraction is correct") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_sub(c, a, b);
			fp2_sub_integ(d, a, b);
			TEST_ASSERT(fp2_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	fp2_free(d);
	return code;
}

static int doubling2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c;

	fp2_null(a);
	fp2_null(b);
	fp2_null(c);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);

		TEST_BEGIN("doubling is correct") {
			fp2_rand(a);
			fp2_dbl(b, a);
			fp2_add(c, a, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_QDR == BASIC || !defined(STRIP)
		TEST_BEGIN("basic doubling is correct") {
			fp2_rand(a);
			fp2_dbl(b, a);
			fp2_dbl_basic(c, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_QDR == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated doubling is correct") {
			fp2_rand(a);
			fp2_dbl(b, a);
			fp2_dbl_integ(c, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	return code;
}

static int multiplication2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c, d, e, f;
	bn_t g;

	bn_null(g);
	fp2_null(a);
	fp2_null(b);
	fp2_null(c);
	fp2_null(d);
	fp2_null(e);
	fp2_null(f);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);
		fp2_new(d);
		fp2_new(e);
		fp2_new(f);
		bn_new(g);

		TEST_BEGIN("multiplication is commutative") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_mul(d, a, b);
			fp2_mul(e, b, a);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is associative") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_rand(c);
			fp2_mul(d, a, b);
			fp2_mul(d, d, c);
			fp2_mul(e, b, c);
			fp2_mul(e, a, e);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is distributive") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_rand(c);
			fp2_add(d, a, b);
			fp2_mul(d, c, d);
			fp2_mul(e, c, a);
			fp2_mul(f, c, b);
			fp2_add(e, e, f);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has identity") {
			fp2_rand(a);
			fp2_zero(d);
			fp_set_dig(d[0], 1);
			fp2_mul(e, a, d);
			TEST_ASSERT(fp2_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has zero property") {
			fp2_rand(a);
			fp2_zero(d);
			fp2_mul(e, a, d);
			TEST_ASSERT(fp2_is_zero(e), end);
		} TEST_END;

#if PP_QDR == BASIC || !defined(STRIP)
		TEST_BEGIN("basic multiplication is correct") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_mul(d, a, b);
			fp2_mul_basic(e, b, a);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_QDR == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated multiplication is correct") {
			fp2_rand(a);
			fp2_rand(b);
			fp2_mul(d, a, b);
			fp2_mul_integ(e, b, a);
			TEST_ASSERT(fp2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("multiplication by adjoined root is correct") {
			fp2_rand(a);
			fp2_zero(b);
			fp_set_dig(b[1], 1);
			fp2_mul(c, a, b);
			fp2_mul_art(d, a);
			TEST_ASSERT(fp2_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication by quadratic/cubic non-residue is correct") {
			fp2_rand(a);
			fp2_mul_nor(b, a);
			switch (fp_prime_get_mod8()) {
				case 5:
					fp2_mul_art(c, a);
					break;
				case 3:
					fp_set_dig(c[0], 1);
					fp_set_dig(c[1], 1);
					fp2_mul(c, a, c);
					break;
				case 7:
					fp_set_dig(c[0], fp2_field_get_qnr());
					fp_set_dig(c[1], 1);
					fp2_mul(c, a, c);
					break;
				default:
					fp2_copy(c, b);
					break;
			}
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		}
		TEST_END;

#if PP_QDR == BASIC || !defined(STRIP)
		TEST_BEGIN("basic multiplication by qnr/cnr is correct") {
			fp2_rand(a);
			fp2_mul_nor(b, a);
			fp2_mul_nor_basic(c, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		}
		TEST_END;
#endif

#if PP_QDR == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated multiplication by qnr/cnr is correct") {
			fp2_rand(a);
			fp2_mul_nor(b, a);
			fp2_mul_nor_integ(c, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		}
		TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	fp2_free(d);
	fp2_free(e);
	fp2_free(f);
	bn_free(g);
	return code;
}

static int squaring2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c;

	fp2_null(a);
	fp2_null(b);
	fp2_null(c);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);

		TEST_BEGIN("squaring is correct") {
			fp2_rand(a);
			fp2_mul(b, a, a);
			fp2_sqr(c, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_QDR == BASIC || !defined(STRIP)
		TEST_BEGIN("basic squaring is correct") {
			fp2_rand(a);
			fp2_sqr(b, a);
			fp2_sqr_basic(c, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_QDR == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated squaring is correct") {
			fp2_rand(a);
			fp2_sqr(b, a);
			fp2_sqr_integ(c, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	return code;
}

static int inversion2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c, d[2];

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);
		fp2_new(d[0]);
		fp2_new(d[1]);

		TEST_BEGIN("inversion is correct") {
			do {
				fp2_rand(a);
			} while (fp2_is_zero(a));
			fp2_inv(b, a);
			fp2_mul(c, a, b);
			fp2_zero(b);
			fp_set_dig(b[0], 1);
			TEST_ASSERT(fp2_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("inversion of a unitary element is correct") {
			do {
				fp2_rand(a);
			} while (fp2_is_zero(a));
			fp2_conv_uni(a, a);
			fp2_inv(b, a);
			fp2_inv_uni(c, a);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("simultaneous inversion is correct") {
			do {
				fp2_rand(a);
				fp2_rand(b);
			} while (fp2_is_zero(a) || fp2_is_zero(b));
			fp2_copy(d[0], a);
			fp2_copy(d[1], b);
			fp2_inv(a, a);
			fp2_inv(b, b);
			fp2_inv_sim(d, d, 2);
			TEST_ASSERT(fp2_cmp(d[0], a) == RLC_EQ &&
					fp2_cmp(d[1], b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	fp2_free(d[0]);
	fp2_free(d[1]);
	return code;
}

static int exponentiation2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c;
	bn_t d;

	fp2_null(a);
	fp2_null(b);
	fp2_null(c);
	bn_null(d);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);
		bn_new(d);

		TEST_BEGIN("exponentiation is correct") {
			fp2_rand(a);
			bn_zero(d);
			fp2_exp(c, a, d);
			TEST_ASSERT(fp2_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(d, 1);
			fp2_exp(c, a, d);
			TEST_ASSERT(fp2_cmp(c, a) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_FP_BITS);
			fp2_exp(b, a, d);
			bn_neg(d, d);
			fp2_exp(c, a, d);
			fp2_inv(c, c);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_DIG);
			fp2_exp_dig(b, a, d->dp[0]);
			fp2_exp(c, a, d);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and exponentiation are consistent") {
			fp2_rand(a);
			fp2_frb(b, a, 0);
			TEST_ASSERT(fp2_cmp(a, b) == RLC_EQ, end);
			fp2_frb(b, a, 1);
			d->sign = RLC_POS;
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			fp2_exp(c, a, d);
			TEST_ASSERT(fp2_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and squared frobenius are consistent") {
			fp2_rand(a);
			fp2_frb(b, a, 1);
			fp2_frb(b, b, 1);
			fp2_frb(c, a, 2);
			TEST_ASSERT(fp2_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

        TEST_BEGIN("exponentiation of unitary element is correct") {
			bn_rand(d, RLC_POS, RLC_FP_BITS);
			fp2_rand(a);
			fp2_conv_uni(a, a);
			fp2_exp(b, a, d);
			fp2_exp_uni(c, a, d);
			TEST_ASSERT(fp2_cmp(b, c) == RLC_EQ, end);
        } TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	bn_free(d);
	return code;
}

static int compression2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c, d;
	uint8_t bin[2 * RLC_FP_BYTES];

	fp2_null(a);
	fp2_null(b);
	fp2_null(c);
	fp2_null(d);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);
		fp2_new(d);

		TEST_BEGIN("compression is consistent") {
			fp2_rand(a);
			fp2_pck(b, a);
			TEST_ASSERT(fp2_upk(c, b) == 1, end);
			TEST_ASSERT(fp2_cmp(a, c) == RLC_EQ, end);
			fp2_rand(a);
			fp2_conv_uni(b, a);
			fp2_pck(c, b);
			TEST_ASSERT(fp2_upk(d, c) == 1, end);
			TEST_ASSERT(fp2_cmp(b, d) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("compression is consitent with reading and writing") {
			fp2_rand(a);
			fp2_conv_uni(a, a);
			fp2_write_bin(bin, RLC_FP_BYTES + 1, a, 1);
			fp2_read_bin(b, bin, RLC_FP_BYTES + 1);
			TEST_ASSERT(fp2_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("getting the size of a compressed field element is correct") {
			fp2_rand(a);
			fp2_conv_uni(a, a);
			TEST_ASSERT(fp2_size_bin(a, 0) == 2 * RLC_FP_BYTES, end);
			TEST_ASSERT(fp2_size_bin(a, 1) == RLC_FP_BYTES + 1, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	fp2_free(d);
	return code;
}

static int square_root2(void) {
	int code = RLC_ERR;
	fp2_t a, b, c;
	int r;

	fp2_null(a);
	fp2_null(b);
	fp2_null(c);

	TRY {
		fp2_new(a);
		fp2_new(b);
		fp2_new(c);

		TEST_BEGIN("square root extraction is correct") {
			fp2_rand(a);
			fp2_sqr(c, a);
			r = fp2_srt(b, c);
			fp2_neg(c, b);
			TEST_ASSERT(r, end);
			TEST_ASSERT(fp2_cmp(b, a) == RLC_EQ ||
					fp2_cmp(c, a) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp2_free(a);
	fp2_free(b);
	fp2_free(c);
	return code;
}

static int memory3(void) {
	err_t e;
	int code = RLC_ERR;
	fp3_t a;

	fp3_null(a);

	TRY {
		TEST_BEGIN("memory can be allocated") {
			fp3_new(a);
			fp3_free(a);
		} TEST_END;
	} CATCH(e) {
		switch (e) {
			case ERR_NO_MEMORY:
				util_print("FATAL ERROR!\n");
				ERROR(end);
				break;
		}
	}
	(void)a;
	code = RLC_OK;
  end:
	return code;
}

static int util3(void) {
	int code = RLC_ERR;
	uint8_t bin[3 * RLC_FP_BYTES];
	fp3_t a, b, c;
	dig_t d;

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);

		TEST_BEGIN("comparison is consistent") {
			fp3_rand(a);
			fp3_rand(b);
			if (fp3_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp3_cmp(b, a) == RLC_NE, end);
			}
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_rand(c);
			if (fp3_cmp(a, c) != RLC_EQ) {
				fp3_copy(c, a);
				TEST_ASSERT(fp3_cmp(c, a) == RLC_EQ, end);
			}
			if (fp3_cmp(b, c) != RLC_EQ) {
				fp3_copy(c, b);
				TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
			}
		}
		TEST_END;

		TEST_BEGIN("negation is consistent") {
			fp3_rand(a);
			fp3_neg(b, a);
			if (fp3_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp3_cmp(b, a) == RLC_NE, end);
			}
			fp3_neg(b, b);
			TEST_ASSERT(fp3_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and comparison are consistent") {
			do {
				fp3_rand(a);
			} while (fp3_is_zero(a));
			fp3_zero(c);
			TEST_ASSERT(fp3_cmp(a, c) == RLC_NE, end);
			TEST_ASSERT(fp3_cmp(c, a) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to random and comparison are consistent") {
			do {
				fp3_rand(a);
			} while (fp3_is_zero(a));
			fp3_zero(c);
			TEST_ASSERT(fp3_cmp(a, c) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and zero test are consistent") {
			fp3_zero(a);
			TEST_ASSERT(fp3_is_zero(a), end);
		}
		TEST_END;

		TEST_BEGIN("assignment to a constant and comparison are consistent") {
			rand_bytes((uint8_t *)&d, (RLC_DIG / 8));
			fp3_set_dig(a, d);
			TEST_ASSERT(fp3_cmp_dig(a, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("reading and writing a finite field element are consistent") {
			fp3_rand(a);
			fp3_write_bin(bin, sizeof(bin), a);
			fp3_read_bin(b, bin, sizeof(bin));
			TEST_ASSERT(fp3_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("getting the size of a prime field element is correct") {
			fp3_rand(a);
			TEST_ASSERT(fp3_size_bin(a) == 3 * RLC_FP_BYTES, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	return code;
}

static int addition3(void) {
	int code = RLC_ERR;
	fp3_t a, b, c, d, e;

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);
	fp3_null(d);
	fp3_null(e);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);
		fp3_new(d);
		fp3_new(e);

		TEST_BEGIN("addition is commutative") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_add(d, a, b);
			fp3_add(e, b, a);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition is associative") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_rand(c);
			fp3_add(d, a, b);
			fp3_add(d, d, c);
			fp3_add(e, b, c);
			fp3_add(e, a, e);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has identity") {
			fp3_rand(a);
			fp3_zero(d);
			fp3_add(e, a, d);
			TEST_ASSERT(fp3_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has inverse") {
			fp3_rand(a);
			fp3_neg(d, a);
			fp3_add(e, a, d);
			TEST_ASSERT(fp3_is_zero(e), end);
		} TEST_END;

#if PP_CBC == BASIC || !defined(STRIP)
		TEST_BEGIN("basic addition is correct") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_add(d, a, b);
			fp3_add_basic(e, a, b);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_CBC == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated addition is correct") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_add(d, a, b);
			fp3_add_integ(e, a, b);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	fp3_free(d);
	fp3_free(e);
	return code;
}

static int subtraction3(void) {
	int code = RLC_ERR;
	fp3_t a, b, c, d;

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);
	fp3_null(d);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);
		fp3_new(d);

		TEST_BEGIN("subtraction is anti-commutative") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_sub(c, a, b);
			fp3_sub(d, b, a);
			fp3_neg(d, d);
			TEST_ASSERT(fp3_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has identity") {
			fp3_rand(a);
			fp3_zero(c);
			fp3_sub(d, a, c);
			TEST_ASSERT(fp3_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has inverse") {
			fp3_rand(a);
			fp3_sub(c, a, a);
			TEST_ASSERT(fp3_is_zero(c), end);
		}
		TEST_END;

#if PP_CBC == BASIC || !defined(STRIP)
		TEST_BEGIN("basic subtraction is correct") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_sub(c, a, b);
			fp3_sub_basic(d, a, b);
			TEST_ASSERT(fp3_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_CBC == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated subtraction is correct") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_sub(c, a, b);
			fp3_sub_integ(d, a, b);
			TEST_ASSERT(fp3_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	fp3_free(d);
	return code;
}

static int doubling3(void) {
	int code = RLC_ERR;
	fp3_t a, b, c;

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);

		TEST_BEGIN("doubling is correct") {
			fp3_rand(a);
			fp3_dbl(b, a);
			fp3_add(c, a, a);
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_CBC == BASIC || !defined(STRIP)
		TEST_BEGIN("basic doubling is correct") {
			fp3_rand(a);
			fp3_dbl(b, a);
			fp3_dbl_basic(c, a);
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_CBC == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated doubling is correct") {
			fp3_rand(a);
			fp3_dbl(b, a);
			fp3_dbl_integ(c, a);
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	return code;
}

static int multiplication3(void) {
	int code = RLC_ERR;
	fp3_t a, b, c, d, e, f;
	bn_t g;

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);
	fp3_null(d);
	fp3_null(e);
	fp3_null(f);
	bn_null(g);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);
		fp3_new(d);
		fp3_new(e);
		fp3_new(f);
		bn_new(g);

		TEST_BEGIN("multiplication is commutative") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_mul(d, a, b);
			fp3_mul(e, b, a);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is associative") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_rand(c);
			fp3_mul(d, a, b);
			fp3_mul(d, d, c);
			fp3_mul(e, b, c);
			fp3_mul(e, a, e);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is distributive") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_rand(c);
			fp3_add(d, a, b);
			fp3_mul(d, c, d);
			fp3_mul(e, c, a);
			fp3_mul(f, c, b);
			fp3_add(e, e, f);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has identity") {
			fp3_rand(a);
			fp3_zero(d);
			fp_set_dig(d[0], 1);
			fp3_mul(e, a, d);
			TEST_ASSERT(fp3_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has zero property") {
			fp3_rand(a);
			fp3_zero(d);
			fp3_mul(e, a, d);
			TEST_ASSERT(fp3_is_zero(e), end);
		} TEST_END;

		TEST_BEGIN("multiplication by adjoined root is correct") {
			fp3_rand(a);
			fp3_zero(b);
			fp_set_dig(b[1], 1);
			fp3_mul(c, a, b);
			fp3_mul_art(d, a);
			TEST_ASSERT(fp3_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication by quadratic/cubic non-residue is correct") {
			fp3_rand(a);
			fp3_mul_nor(b, a);
			fp3_mul_art(c, a);
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		}
		TEST_END;

#if PP_CBC == BASIC || !defined(STRIP)
		TEST_BEGIN("basic multiplication is correct") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_mul(d, a, b);
			fp3_mul_basic(e, b, a);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_CBC == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated multiplication is correct") {
			fp3_rand(a);
			fp3_rand(b);
			fp3_mul(d, a, b);
			fp3_mul_integ(e, b, a);
			TEST_ASSERT(fp3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("multiplication by adjoined root is correct") {
			fp3_rand(a);
			fp3_zero(b);
			fp_set_dig(b[1], 1);
			fp3_mul(c, a, b);
			fp3_mul_art(d, a);
			TEST_ASSERT(fp3_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication by quadratic/cubic non-residue is correct") {
			fp3_rand(a);
			fp3_mul_nor(b, a);
			switch (fp_prime_get_mod8()) {
				case 5:
					fp3_mul_art(c, a);
					break;
				case 3:
					fp3_mul_art(c, a);
					fp3_add(c, c, a);
					break;
				case 7:
					fp3_mul_art(d, a);
					fp3_dbl(c, a);
					fp_prime_back(g, ep_curve_get_b());
					for (int i = 1; i < bn_bits(g) / 2; i++) {
						fp3_dbl(c, c);
					}
					fp3_add(c, c, d);
					break;
			}
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	fp3_free(d);
	fp3_free(e);
	fp3_free(f);
	bn_free(g);
	return code;
}

static int squaring3(void) {
	int code = RLC_ERR;
	fp3_t a, b, c;

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);

		TEST_BEGIN("squaring is correct") {
			fp3_rand(a);
			fp3_mul(b, a, a);
			fp3_sqr(c, a);
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_CBC == BASIC || !defined(STRIP)
		TEST_BEGIN("basic squaring is correct") {
			fp3_rand(a);
			fp3_sqr(b, a);
			fp3_sqr_basic(c, a);
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_CBC == INTEG || !defined(STRIP)
		TEST_BEGIN("integrated squaring is correct") {
			fp3_rand(a);
			fp3_sqr(b, a);
			fp3_sqr_integ(c, a);
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	return code;
}

static int inversion3(void) {
	int code = RLC_ERR;
	fp3_t a, b, c, d[2];

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);
	fp3_null(d[0]);
	fp3_null(d[1]);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);
		fp3_new(d[0]);
		fp3_new(d[1]);

		TEST_BEGIN("inversion is correct") {
			do {
				fp3_rand(a);
			} while (fp3_is_zero(a));
			fp3_inv(b, a);
			fp3_mul(c, a, b);
			fp3_zero(b);
			fp_set_dig(b[0], 1);
			TEST_ASSERT(fp3_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("simultaneous inversion is correct") {
			do {
				fp3_rand(a);
				fp3_rand(b);
			} while (fp2_is_zero(a) || fp2_is_zero(b));
			fp3_copy(d[0], a);
			fp3_copy(d[1], b);
			fp3_inv(a, a);
			fp3_inv(b, b);
			fp3_inv_sim(d, d, 2);
			TEST_ASSERT(fp3_cmp(d[0], a) == RLC_EQ &&
					fp3_cmp(d[1], b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	fp3_free(d[0]);
	fp3_free(d[1]);
	return code;
}

static int exponentiation3(void) {
	int code = RLC_ERR;
	fp3_t a, b, c;
	bn_t d;

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);
	bn_null(d);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);
		bn_new(d);

		TEST_BEGIN("exponentiation is correct") {
			fp3_rand(a);
			bn_zero(d);
			fp3_exp(c, a, d);
			TEST_ASSERT(fp3_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(d, 1);
			fp3_exp(c, a, d);
			TEST_ASSERT(fp3_cmp(c, a) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_FP_BITS);
			fp3_exp(b, a, d);
			bn_neg(d, d);
			fp3_exp(c, a, d);
			fp3_inv(c, c);
			TEST_ASSERT(fp3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and exponentiation are consistent") {
			fp3_rand(a);
			fp3_frb(b, a, 0);
			TEST_ASSERT(fp3_cmp(a, b) == RLC_EQ, end);
			fp3_frb(b, a, 1);
			d->sign = RLC_POS;
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			fp3_exp(c, a, d);
			TEST_ASSERT(fp3_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and squared frobenius are consistent") {
			fp3_rand(a);
			fp3_frb(b, a, 1);
			fp3_frb(b, b, 1);
			fp3_frb(c, a, 2);
			TEST_ASSERT(fp3_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and cubed frobenius are consistent") {
			fp3_rand(a);
			fp3_frb(b, a, 1);
			fp3_frb(b, b, 1);
			fp3_frb(b, b, 1);
			fp3_frb(c, a, 3);
			TEST_ASSERT(fp3_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	bn_free(d);
	return code;
}

static int square_root3(void) {
	int code = RLC_ERR;
	fp3_t a, b, c;
	int r;

	fp3_null(a);
	fp3_null(b);
	fp3_null(c);

	TRY {
		fp3_new(a);
		fp3_new(b);
		fp3_new(c);

		TEST_BEGIN("square root extraction is correct") {
			fp3_rand(a);
			fp3_sqr(c, a);
			r = fp3_srt(b, c);
			fp3_neg(c, b);
			TEST_ASSERT(r == 1, end);
			TEST_ASSERT(fp3_cmp(b, a) == RLC_EQ || fp3_cmp(c, a) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp3_free(a);
	fp3_free(b);
	fp3_free(c);
	return code;
}

static int memory4(void) {
	err_t e;
	int code = RLC_ERR;
	fp4_t a;

	fp4_null(a);

	TRY {
		TEST_BEGIN("memory can be allocated") {
			fp4_new(a);
			fp4_free(a);
		} TEST_END;
	} CATCH(e) {
		switch (e) {
			case ERR_NO_MEMORY:
				util_print("FATAL ERROR!\n");
				ERROR(end);
				break;
		}
	}
	(void)a;
	code = RLC_OK;
  end:
	return code;
}

static int util4(void) {
	int code = RLC_ERR;
	uint8_t bin[4 * RLC_FP_BYTES];
	fp4_t a, b, c;
	dig_t d;

	fp4_null(a);
	fp4_null(b);
	fp4_null(c);

	TRY {
		fp4_new(a);
		fp4_new(b);
		fp4_new(c);

		TEST_BEGIN("comparison is consistent") {
			fp4_rand(a);
			fp4_rand(b);
			if (fp4_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp4_cmp(b, a) == RLC_NE, end);
			}
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_rand(c);
			if (fp4_cmp(a, c) != RLC_EQ) {
				fp4_copy(c, a);
				TEST_ASSERT(fp4_cmp(c, a) == RLC_EQ, end);
			}
			if (fp4_cmp(b, c) != RLC_EQ) {
				fp4_copy(c, b);
				TEST_ASSERT(fp4_cmp(b, c) == RLC_EQ, end);
			}
		}
		TEST_END;

		TEST_BEGIN("negation is consistent") {
			fp4_rand(a);
			fp4_neg(b, a);
			if (fp4_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp4_cmp(b, a) == RLC_NE, end);
			}
			fp4_neg(b, b);
			TEST_ASSERT(fp4_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and comparison are consistent") {
			do {
				fp4_rand(a);
			} while (fp4_is_zero(a));
			fp4_zero(c);
			TEST_ASSERT(fp4_cmp(a, c) == RLC_NE, end);
			TEST_ASSERT(fp4_cmp(c, a) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to random and comparison are consistent") {
			do {
				fp4_rand(a);
			} while (fp4_is_zero(a));
			fp4_zero(c);
			TEST_ASSERT(fp4_cmp(a, c) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and zero test are consistent") {
			fp4_zero(a);
			TEST_ASSERT(fp4_is_zero(a), end);
		}
		TEST_END;

		TEST_BEGIN("assignment to a constant and comparison are consistent") {
			rand_bytes((uint8_t *)&d, (RLC_DIG / 8));
			fp4_set_dig(a, d);
			TEST_ASSERT(fp4_cmp_dig(a, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("reading and writing a finite field element are consistent") {
			fp4_rand(a);
			fp4_write_bin(bin, sizeof(bin), a);
			fp4_read_bin(b, bin, sizeof(bin));
			TEST_ASSERT(fp4_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("getting the size of a prime field element is correct") {
			fp4_rand(a);
			TEST_ASSERT(fp4_size_bin(a) == 4 * RLC_FP_BYTES, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp4_free(a);
	fp4_free(b);
	fp4_free(c);
	return code;
}

static int addition4(void) {
	int code = RLC_ERR;
	fp4_t a, b, c, d, e;

	fp4_null(a);
	fp4_null(b);
	fp4_null(c);
	fp4_null(d);
	fp4_null(e);

	TRY {
		fp4_new(a);
		fp4_new(b);
		fp4_new(c);
		fp4_new(d);
		fp4_new(e);

		TEST_BEGIN("addition is commutative") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_add(d, a, b);
			fp4_add(e, b, a);
			TEST_ASSERT(fp4_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition is associative") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_rand(c);
			fp4_add(d, a, b);
			fp4_add(d, d, c);
			fp4_add(e, b, c);
			fp4_add(e, a, e);
			TEST_ASSERT(fp4_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has identity") {
			fp4_rand(a);
			fp4_zero(d);
			fp4_add(e, a, d);
			TEST_ASSERT(fp4_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has inverse") {
			fp4_rand(a);
			fp4_neg(d, a);
			fp4_add(e, a, d);
			TEST_ASSERT(fp4_is_zero(e), end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp4_free(a);
	fp4_free(b);
	fp4_free(c);
	fp4_free(d);
	fp4_free(e);
	return code;
}

static int subtraction4(void) {
	int code = RLC_ERR;
	fp4_t a, b, c, d;

	fp4_null(a);
	fp4_null(b);
	fp4_null(c);
	fp4_null(d);

	TRY {
		fp4_new(a);
		fp4_new(b);
		fp4_new(c);
		fp4_new(d);

		TEST_BEGIN("subtraction is anti-commutative") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_sub(c, a, b);
			fp4_sub(d, b, a);
			fp4_neg(d, d);
			TEST_ASSERT(fp4_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has identity") {
			fp4_rand(a);
			fp4_zero(c);
			fp4_sub(d, a, c);
			TEST_ASSERT(fp4_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has inverse") {
			fp4_rand(a);
			fp4_sub(c, a, a);
			TEST_ASSERT(fp4_is_zero(c), end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp4_free(a);
	fp4_free(b);
	fp4_free(c);
	fp4_free(d);
	return code;
}

static int doubling4(void) {
	int code = RLC_ERR;
	fp4_t a, b, c;

	fp4_null(a);
	fp4_null(b);
	fp4_null(c);

	TRY {
		fp4_new(a);
		fp4_new(b);
		fp4_new(c);

		TEST_BEGIN("doubling is correct") {
			fp4_rand(a);
			fp4_dbl(b, a);
			fp4_add(c, a, a);
			TEST_ASSERT(fp4_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp4_free(a);
	fp4_free(b);
	fp4_free(c);
	return code;
}

static int multiplication4(void) {
	int code = RLC_ERR;
	fp4_t a, b, c, d, e, f;

	fp4_null(a);
	fp4_null(b);
	fp4_null(c);
	fp4_null(d);
	fp4_null(e);
	fp4_null(f);

	TRY {
		fp4_new(a);
		fp4_new(b);
		fp4_new(c);
		fp4_new(d);
		fp4_new(e);
		fp4_new(f);

		TEST_BEGIN("multiplication is commutative") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_mul(d, a, b);
			fp4_mul(e, b, a);
			TEST_ASSERT(fp4_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is associative") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_rand(c);
			fp4_mul(d, a, b);
			fp4_mul(d, d, c);
			fp4_mul(e, b, c);
			fp4_mul(e, a, e);
			TEST_ASSERT(fp4_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is distributive") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_rand(c);
			fp4_add(d, a, b);
			fp4_mul(d, c, d);
			fp4_mul(e, c, a);
			fp4_mul(f, c, b);
			fp4_add(e, e, f);
			TEST_ASSERT(fp4_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has identity") {
			fp4_zero(d);
			fp_set_dig(d[0][0], 1);
			fp4_mul(e, a, d);
			TEST_ASSERT(fp4_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has zero property") {
			fp4_zero(d);
			fp4_mul(e, a, d);
			TEST_ASSERT(fp4_is_zero(e), end);
		} TEST_END;

		TEST_BEGIN("multiplication by adjoined root is correct") {
			fp4_rand(a);
			fp4_zero(b);
			fp_set_dig(b[1][0], 1);
			fp4_mul(c, a, b);
			fp4_mul_art(d, a);
			TEST_ASSERT(fp4_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic multiplication is correct") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_mul(c, a, b);
			fp4_mul_basic(d, a, b);
			TEST_ASSERT(fp4_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced multiplication is correct") {
			fp4_rand(a);
			fp4_rand(b);
			fp4_mul(c, a, b);
			fp4_mul_lazyr(d, a, b);
			TEST_ASSERT(fp4_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp4_free(a);
	fp4_free(b);
	fp4_free(c);
	fp4_free(d);
	fp4_free(e);
	fp4_free(f);
	return code;
}

static int squaring4(void) {
	int code = RLC_ERR;
	fp4_t a, b, c;

	fp4_null(a);
	fp4_null(b);
	fp4_null(c);

	TRY {
		fp4_new(a);
		fp4_new(b);
		fp4_new(c);

		TEST_BEGIN("squaring is correct") {
			fp4_rand(a);
			fp4_mul(b, a, a);
			fp4_sqr(c, a);
			TEST_ASSERT(fp4_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic squaring is correct") {
			fp4_rand(a);
			fp4_sqr(b, a);
			fp4_sqr_basic(c, a);
			TEST_ASSERT(fp4_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced squaring is correct") {
			fp4_rand(a);
			fp4_sqr(b, a);
			fp4_sqr_lazyr(c, a);
			TEST_ASSERT(fp4_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp4_free(a);
	fp4_free(b);
	fp4_free(c);
	return code;
}

static int inversion4(void) {
	int code = RLC_ERR;
	fp4_t a, b, c;

	fp4_null(a);
	fp4_null(b);
	fp4_null(c);

	TRY {
		fp4_new(a);
		fp4_new(b);
		fp4_new(c);

		TEST_BEGIN("inversion is correct") {
			do {
				fp4_rand(a);
			} while (fp4_is_zero(a));
			fp4_inv(b, a);
			fp4_mul(c, a, b);
			fp4_zero(b);
			fp_set_dig(b[0][0], 1);
			TEST_ASSERT(fp4_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp4_free(a);
	fp4_free(b);
	fp4_free(c);
	return code;
}

static int exponentiation4(void) {
	int code = RLC_ERR;
	fp4_t a, b, c;
	bn_t d;

	fp4_null(a);
	fp4_null(b);
	fp4_null(c);
	bn_null(d);

	TRY {
		fp4_new(a);
		fp4_new(b);
		fp4_new(c);
		bn_new(d);

		TEST_BEGIN("exponentiation is correct") {
			fp4_rand(a);
			bn_zero(d);
			fp4_exp(c, a, d);
			TEST_ASSERT(fp4_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(d, 1);
			fp4_exp(c, a, d);
			TEST_ASSERT(fp4_cmp(c, a) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_FP_BITS);
			fp4_exp(b, a, d);
			bn_neg(d, d);
			fp4_exp(c, a, d);
			fp4_inv(c, c);
			TEST_ASSERT(fp4_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and exponentiation are consistent") {
			fp4_rand(a);
			fp4_frb(b, a, 0);
			TEST_ASSERT(fp4_cmp(a, b) == RLC_EQ, end);
			fp4_frb(b, a, 1);
			d->sign = RLC_POS;
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			fp4_exp(c, a, d);
			TEST_ASSERT(fp4_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and squared frobenius are consistent") {
			fp4_rand(a);
			fp4_frb(b, a, 1);
			fp4_frb(b, b, 1);
			fp4_frb(c, a, 2);
			TEST_ASSERT(fp4_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp4_free(a);
	fp4_free(b);
	fp4_free(c);
	bn_free(d);
	return code;
}

static int memory6(void) {
	err_t e;
	int code = RLC_ERR;
	fp6_t a;

	fp6_null(a);

	TRY {
		TEST_BEGIN("memory can be allocated") {
			fp6_new(a);
			fp6_free(a);
		} TEST_END;
	} CATCH(e) {
		switch (e) {
			case ERR_NO_MEMORY:
				util_print("FATAL ERROR!\n");
				ERROR(end);
				break;
		}
	}
	(void)a;
	code = RLC_OK;
  end:
	return code;
}

static int util6(void) {
	int code = RLC_ERR;
	uint8_t bin[6 * RLC_FP_BYTES];
	fp6_t a, b, c;
	dig_t d;

	fp6_null(a);
	fp6_null(b);
	fp6_null(c);

	TRY {
		fp6_new(a);
		fp6_new(b);
		fp6_new(c);

		TEST_BEGIN("comparison is consistent") {
			fp6_rand(a);
			fp6_rand(b);
			if (fp6_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp6_cmp(b, a) == RLC_NE, end);
			}
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_rand(c);
			if (fp6_cmp(a, c) != RLC_EQ) {
				fp6_copy(c, a);
				TEST_ASSERT(fp6_cmp(c, a) == RLC_EQ, end);
			}
			if (fp6_cmp(b, c) != RLC_EQ) {
				fp6_copy(c, b);
				TEST_ASSERT(fp6_cmp(b, c) == RLC_EQ, end);
			}
		}
		TEST_END;

		TEST_BEGIN("negation is consistent") {
			fp6_rand(a);
			fp6_neg(b, a);
			if (fp6_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp6_cmp(b, a) == RLC_NE, end);
			}
			fp6_neg(b, b);
			TEST_ASSERT(fp6_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and comparison are consistent") {
			do {
				fp6_rand(a);
			} while (fp6_is_zero(a));
			fp6_zero(c);
			TEST_ASSERT(fp6_cmp(a, c) == RLC_NE, end);
			TEST_ASSERT(fp6_cmp(c, a) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to random and comparison are consistent") {
			do {
				fp6_rand(a);
			} while (fp6_is_zero(a));
			fp6_zero(c);
			TEST_ASSERT(fp6_cmp(a, c) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and zero test are consistent") {
			fp6_zero(a);
			TEST_ASSERT(fp6_is_zero(a), end);
		}
		TEST_END;

		TEST_BEGIN("assignment to a constant and comparison are consistent") {
			rand_bytes((uint8_t *)&d, (RLC_DIG / 8));
			fp6_set_dig(a, d);
			TEST_ASSERT(fp6_cmp_dig(a, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("reading and writing a finite field element are consistent") {
			fp6_rand(a);
			fp6_write_bin(bin, sizeof(bin), a);
			fp6_read_bin(b, bin, sizeof(bin));
			TEST_ASSERT(fp6_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("getting the size of a prime field element is correct") {
			fp6_rand(a);
			TEST_ASSERT(fp6_size_bin(a) == 6 * RLC_FP_BYTES, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp6_free(a);
	fp6_free(b);
	fp6_free(c);
	return code;
}

static int addition6(void) {
	int code = RLC_ERR;
	fp6_t a, b, c, d, e;

	fp6_null(a);
	fp6_null(b);
	fp6_null(c);
	fp6_null(d);
	fp6_null(e);

	TRY {
		fp6_new(a);
		fp6_new(b);
		fp6_new(c);
		fp6_new(d);
		fp6_new(e);

		TEST_BEGIN("addition is commutative") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_add(d, a, b);
			fp6_add(e, b, a);
			TEST_ASSERT(fp6_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition is associative") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_rand(c);
			fp6_add(d, a, b);
			fp6_add(d, d, c);
			fp6_add(e, b, c);
			fp6_add(e, a, e);
			TEST_ASSERT(fp6_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has identity") {
			fp6_rand(a);
			fp6_zero(d);
			fp6_add(e, a, d);
			TEST_ASSERT(fp6_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has inverse") {
			fp6_rand(a);
			fp6_neg(d, a);
			fp6_add(e, a, d);
			TEST_ASSERT(fp6_is_zero(e), end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp6_free(a);
	fp6_free(b);
	fp6_free(c);
	fp6_free(d);
	fp6_free(e);
	return code;
}

static int subtraction6(void) {
	int code = RLC_ERR;
	fp6_t a, b, c, d;

	fp6_null(a);
	fp6_null(b);
	fp6_null(c);
	fp6_null(d);

	TRY {
		fp6_new(a);
		fp6_new(b);
		fp6_new(c);
		fp6_new(d);

		TEST_BEGIN("subtraction is anti-commutative") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_sub(c, a, b);
			fp6_sub(d, b, a);
			fp6_neg(d, d);
			TEST_ASSERT(fp6_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has identity") {
			fp6_rand(a);
			fp6_zero(c);
			fp6_sub(d, a, c);
			TEST_ASSERT(fp6_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has inverse") {
			fp6_rand(a);
			fp6_sub(c, a, a);
			TEST_ASSERT(fp6_is_zero(c), end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp6_free(a);
	fp6_free(b);
	fp6_free(c);
	fp6_free(d);
	return code;
}

static int doubling6(void) {
	int code = RLC_ERR;
	fp6_t a, b, c;

	fp6_null(a);
	fp6_null(b);
	fp6_null(c);

	TRY {
		fp6_new(a);
		fp6_new(b);
		fp6_new(c);

		TEST_BEGIN("doubling is correct") {
			fp6_rand(a);
			fp6_dbl(b, a);
			fp6_add(c, a, a);
			TEST_ASSERT(fp6_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp6_free(a);
	fp6_free(b);
	fp6_free(c);
	return code;
}

static int multiplication6(void) {
	int code = RLC_ERR;
	fp6_t a, b, c, d, e, f;

	fp6_null(a);
	fp6_null(b);
	fp6_null(c);
	fp6_null(d);
	fp6_null(e);
	fp6_null(f);

	TRY {
		fp6_new(a);
		fp6_new(b);
		fp6_new(c);
		fp6_new(d);
		fp6_new(e);
		fp6_new(f);

		TEST_BEGIN("multiplication is commutative") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_mul(d, a, b);
			fp6_mul(e, b, a);
			TEST_ASSERT(fp6_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is associative") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_rand(c);
			fp6_mul(d, a, b);
			fp6_mul(d, d, c);
			fp6_mul(e, b, c);
			fp6_mul(e, a, e);
			TEST_ASSERT(fp6_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is distributive") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_rand(c);
			fp6_add(d, a, b);
			fp6_mul(d, c, d);
			fp6_mul(e, c, a);
			fp6_mul(f, c, b);
			fp6_add(e, e, f);
			TEST_ASSERT(fp6_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has identity") {
			fp6_zero(d);
			fp_set_dig(d[0][0], 1);
			fp6_mul(e, a, d);
			TEST_ASSERT(fp6_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has zero property") {
			fp6_zero(d);
			fp6_mul(e, a, d);
			TEST_ASSERT(fp6_is_zero(e), end);
		} TEST_END;

		TEST_BEGIN("multiplication by adjoined root is correct") {
			fp6_rand(a);
			fp6_zero(b);
			fp_set_dig(b[1][0], 1);
			fp6_mul(c, a, b);
			fp6_mul_art(d, a);
			TEST_ASSERT(fp6_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic multiplication is correct") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_mul(c, a, b);
			fp6_mul_basic(d, a, b);
			TEST_ASSERT(fp6_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced multiplication is correct") {
			fp6_rand(a);
			fp6_rand(b);
			fp6_mul(c, a, b);
			fp6_mul_lazyr(d, a, b);
			TEST_ASSERT(fp6_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp6_free(a);
	fp6_free(b);
	fp6_free(c);
	fp6_free(d);
	fp6_free(e);
	fp6_free(f);
	return code;
}

static int squaring6(void) {
	int code = RLC_ERR;
	fp6_t a, b, c;

	fp6_null(a);
	fp6_null(b);
	fp6_null(c);

	TRY {
		fp6_new(a);
		fp6_new(b);
		fp6_new(c);

		TEST_BEGIN("squaring is correct") {
			fp6_rand(a);
			fp6_mul(b, a, a);
			fp6_sqr(c, a);
			TEST_ASSERT(fp6_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic squaring is correct") {
			fp6_rand(a);
			fp6_sqr(b, a);
			fp6_sqr_basic(c, a);
			TEST_ASSERT(fp6_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced squaring is correct") {
			fp6_rand(a);
			fp6_sqr(b, a);
			fp6_sqr_lazyr(c, a);
			TEST_ASSERT(fp6_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp6_free(a);
	fp6_free(b);
	fp6_free(c);
	return code;
}

static int inversion6(void) {
	int code = RLC_ERR;
	fp6_t a, b, c;

	fp6_null(a);
	fp6_null(b);
	fp6_null(c);

	TRY {
		fp6_new(a);
		fp6_new(b);
		fp6_new(c);

		TEST_BEGIN("inversion is correct") {
			do {
				fp6_rand(a);
			} while (fp6_is_zero(a));
			fp6_inv(b, a);
			fp6_mul(c, a, b);
			fp6_zero(b);
			fp_set_dig(b[0][0], 1);
			TEST_ASSERT(fp6_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp6_free(a);
	fp6_free(b);
	fp6_free(c);
	return code;
}

static int exponentiation6(void) {
	int code = RLC_ERR;
	fp6_t a, b, c;
	bn_t d;

	fp6_null(a);
	fp6_null(b);
	fp6_null(c);
	bn_null(d);

	TRY {
		fp6_new(a);
		fp6_new(b);
		fp6_new(c);
		bn_new(d);

		TEST_BEGIN("exponentiation is correct") {
			fp6_rand(a);
			bn_zero(d);
			fp6_exp(c, a, d);
			TEST_ASSERT(fp6_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(d, 1);
			fp6_exp(c, a, d);
			TEST_ASSERT(fp6_cmp(c, a) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_FP_BITS);
			fp6_exp(b, a, d);
			bn_neg(d, d);
			fp6_exp(c, a, d);
			fp6_inv(c, c);
			TEST_ASSERT(fp6_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and exponentiation are consistent") {
			fp6_rand(a);
			fp6_frb(b, a, 0);
			TEST_ASSERT(fp6_cmp(a, b) == RLC_EQ, end);
			fp6_frb(b, a, 1);
			d->sign = RLC_POS;
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			fp6_exp(c, a, d);
			TEST_ASSERT(fp6_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and squared frobenius are consistent") {
			fp6_rand(a);
			fp6_frb(b, a, 1);
			fp6_frb(b, b, 1);
			fp6_frb(c, a, 2);
			TEST_ASSERT(fp6_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp6_free(a);
	fp6_free(b);
	fp6_free(c);
	bn_free(d);
	return code;
}

static int memory8(void) {
	err_t e;
	int code = RLC_ERR;
	fp8_t a;

	fp8_null(a);

	TRY {
		TEST_BEGIN("memory can be allocated") {
			fp8_new(a);
			fp8_free(a);
		} TEST_END;
	} CATCH(e) {
		switch (e) {
			case ERR_NO_MEMORY:
				util_print("FATAL ERROR!\n");
				ERROR(end);
				break;
		}
	}
	(void)a;
	code = RLC_OK;
  end:
	return code;
}

static int util8(void) {
	int code = RLC_ERR;
	uint8_t bin[8 * RLC_FP_BYTES];
	fp8_t a, b, c;
	dig_t d;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);

		TEST_BEGIN("comparison is consistent") {
			fp8_rand(a);
			fp8_rand(b);
			if (fp8_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp8_cmp(b, a) == RLC_NE, end);
			}
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_rand(c);
			if (fp8_cmp(a, c) != RLC_EQ) {
				fp8_copy(c, a);
				TEST_ASSERT(fp8_cmp(c, a) == RLC_EQ, end);
			}
			if (fp8_cmp(b, c) != RLC_EQ) {
				fp8_copy(c, b);
				TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
			}
		}
		TEST_END;

		TEST_BEGIN("negation is consistent") {
			fp8_rand(a);
			fp8_neg(b, a);
			if (fp8_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp8_cmp(b, a) == RLC_NE, end);
			}
			fp8_neg(b, b);
			TEST_ASSERT(fp8_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and comparison are consistent") {
			do {
				fp8_rand(a);
			} while (fp8_is_zero(a));
			fp8_zero(c);
			TEST_ASSERT(fp8_cmp(a, c) == RLC_NE, end);
			TEST_ASSERT(fp8_cmp(c, a) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to random and comparison are consistent") {
			do {
				fp8_rand(a);
			} while (fp8_is_zero(a));
			fp8_zero(c);
			TEST_ASSERT(fp8_cmp(a, c) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and zero test are consistent") {
			fp8_zero(a);
			TEST_ASSERT(fp8_is_zero(a), end);
		}
		TEST_END;

		TEST_BEGIN("assignment to a constant and comparison are consistent") {
			rand_bytes((uint8_t *)&d, (RLC_DIG / 8));
			fp8_set_dig(a, d);
			TEST_ASSERT(fp8_cmp_dig(a, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("reading and writing a finite field element are consistent") {
			fp8_rand(a);
			fp8_write_bin(bin, sizeof(bin), a);
			fp8_read_bin(b, bin, sizeof(bin));
			TEST_ASSERT(fp8_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("getting the size of a finite field element is correct") {
			fp8_rand(a);
			TEST_ASSERT(fp8_size_bin(a, 0) == 8 * RLC_FP_BYTES, end);
			fp8_conv_uni(a, a);
			TEST_ASSERT(fp8_size_bin(a, 1) == 4 * RLC_FP_BYTES, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	return code;
}

static int addition8(void) {
	int code = RLC_ERR;
	fp8_t a, b, c, d, e;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);
	fp8_null(d);
	fp8_null(e);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);
		fp8_new(d);
		fp8_new(e);

		TEST_BEGIN("addition is commutative") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_add(d, a, b);
			fp8_add(e, b, a);
			TEST_ASSERT(fp8_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition is associative") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_rand(c);
			fp8_add(d, a, b);
			fp8_add(d, d, c);
			fp8_add(e, b, c);
			fp8_add(e, a, e);
			TEST_ASSERT(fp8_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has identity") {
			fp8_rand(a);
			fp8_zero(d);
			fp8_add(e, a, d);
			TEST_ASSERT(fp8_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has inverse") {
			fp8_rand(a);
			fp8_neg(d, a);
			fp8_add(e, a, d);
			TEST_ASSERT(fp8_is_zero(e), end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	fp8_free(d);
	fp8_free(e);
	return code;
}

static int subtraction8(void) {
	int code = RLC_ERR;
	fp8_t a, b, c, d;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);
	fp8_null(d);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);
		fp8_new(d);

		TEST_BEGIN("subtraction is anti-commutative") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_sub(c, a, b);
			fp8_sub(d, b, a);
			fp8_neg(d, d);
			TEST_ASSERT(fp8_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has identity") {
			fp8_rand(a);
			fp8_zero(c);
			fp8_sub(d, a, c);
			TEST_ASSERT(fp8_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has inverse") {
			fp8_rand(a);
			fp8_sub(c, a, a);
			TEST_ASSERT(fp8_is_zero(c), end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	fp8_free(d);
	return code;
}

static int doubling8(void) {
	int code = RLC_ERR;
	fp8_t a, b, c;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);

		TEST_BEGIN("doubling is correct") {
			fp8_rand(a);
			fp8_dbl(b, a);
			fp8_add(c, a, a);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	return code;
}

static int multiplication8(void) {
	int code = RLC_ERR;
	fp8_t a, b, c, d, e, f;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);
	fp8_null(d);
	fp8_null(e);
	fp8_null(f);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);
		fp8_new(d);
		fp8_new(e);
		fp8_new(f);

		TEST_BEGIN("multiplication is commutative") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_mul(d, a, b);
			fp8_mul(e, b, a);
			TEST_ASSERT(fp8_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is associative") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_rand(c);
			fp8_mul(d, a, b);
			fp8_mul(d, d, c);
			fp8_mul(e, b, c);
			fp8_mul(e, a, e);
			TEST_ASSERT(fp8_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is distributive") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_rand(c);
			fp8_add(d, a, b);
			fp8_mul(d, c, d);
			fp8_mul(e, c, a);
			fp8_mul(f, c, b);
			fp8_add(e, e, f);
			TEST_ASSERT(fp8_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has identity") {
			fp8_zero(d);
			fp_set_dig(d[0][0][0], 1);
			fp8_mul(e, a, d);
			TEST_ASSERT(fp8_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has zero property") {
			fp8_zero(d);
			fp8_mul(e, a, d);
			TEST_ASSERT(fp8_is_zero(e), end);
		} TEST_END;

		TEST_BEGIN("multiplication by adjoined root is correct") {
			fp8_rand(a);
			fp8_zero(b);
			fp_set_dig(b[1][0][0], 1);
			fp8_mul(c, a, b);
			fp8_mul_art(d, a);
			TEST_ASSERT(fp8_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic multiplication is correct") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_mul(c, a, b);
			fp8_mul_basic(d, a, b);
			TEST_ASSERT(fp8_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced multiplication is correct") {
			fp8_rand(a);
			fp8_rand(b);
			fp8_mul(c, a, b);
			fp8_mul_lazyr(d, a, b);
			TEST_ASSERT(fp8_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	fp8_free(d);
	fp8_free(e);
	fp8_free(f);
	return code;
}

static int squaring8(void) {
	int code = RLC_ERR;
	fp8_t a, b, c;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);

		TEST_BEGIN("squaring is correct") {
			fp8_rand(a);
			fp8_mul(b, a, a);
			fp8_sqr(c, a);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic squaring is correct") {
			fp8_rand(a);
			fp8_sqr(b, a);
			fp8_sqr_basic(c, a);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced squaring is correct") {
			fp8_rand(a);
			fp8_sqr(b, a);
			fp8_sqr_lazyr(c, a);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	return code;
}

static int cyclotomic8(void) {
	int code = RLC_ERR;
	fp8_t a, b, c;
	bn_t f;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);
	bn_null(f);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);
		bn_new(f);

		TEST_BEGIN("cyclotomic test is correct") {
			fp8_rand(a);
			fp8_conv_uni(a, a);
			TEST_ASSERT(fp8_test_uni(a) == 1, end);
		} TEST_END;

		TEST_BEGIN("cyclotomic squaring is correct") {
			fp8_rand(a);
			fp8_conv_uni(a, a);
			fp8_sqr(b, a);
			fp8_sqr_uni(c, a);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

        TEST_BEGIN("cyclotomic exponentiation is correct") {
			fp8_rand(a);
			fp8_conv_uni(a, a);
			bn_zero(f);
			fp8_exp_uni(c, a, f);
			TEST_ASSERT(fp8_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(f, 1);
			fp8_exp_uni(c, a, f);
			TEST_ASSERT(fp8_cmp(c, a) == RLC_EQ, end);
			bn_rand(f, RLC_POS, RLC_FP_BITS);
			fp8_exp(b, a, f);
			fp8_exp_uni(c, a, f);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
			bn_rand(f, RLC_POS, RLC_FP_BITS);
			fp8_exp_uni(b, a, f);
			bn_neg(f, f);
			fp8_exp_uni(c, a, f);
			fp8_inv_uni(c, c);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
        } TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	bn_free(f);
	return code;
}

static int inversion8(void) {
	int code = RLC_ERR;
	fp8_t a, b, c;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);

		TEST_BEGIN("inversion is correct") {
			do {
				fp8_rand(a);
			} while (fp8_is_zero(a));
			fp8_inv(b, a);
			fp8_mul(c, a, b);
			fp8_zero(b);
			fp_set_dig(b[0][0][0], 1);
			TEST_ASSERT(fp8_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("inversion of a unitary element is correct") {
			do {
				fp8_rand(a);
			} while (fp8_is_zero(a));
			fp8_conv_uni(a, a);
			fp8_inv(b, a);
			fp8_inv_uni(c, a);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	return code;
}

static int exponentiation8(void) {
	int code = RLC_ERR;
	fp8_t a, b, c;
	bn_t d;

	fp8_null(a);
	fp8_null(b);
	fp8_null(c);
	bn_null(d);

	TRY {
		fp8_new(a);
		fp8_new(b);
		fp8_new(c);
		bn_new(d);

		TEST_BEGIN("exponentiation is correct") {
			fp8_rand(a);
			bn_zero(d);
			fp8_exp(c, a, d);
			TEST_ASSERT(fp8_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(d, 1);
			fp8_exp(c, a, d);
			TEST_ASSERT(fp8_cmp(c, a) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_FP_BITS);
			fp8_exp(b, a, d);
			bn_neg(d, d);
			fp8_exp(c, a, d);
			fp8_inv(c, c);
			TEST_ASSERT(fp8_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		if (fp_prime_get_mod8() == 1) {
			TEST_BEGIN("frobenius and exponentiation are consistent") {
				fp8_rand(a);
				fp8_frb(b, a, 0);
				TEST_ASSERT(fp8_cmp(a, b) == RLC_EQ, end);
				fp8_frb(b, a, 1);
				d->sign = RLC_POS;
				d->used = RLC_FP_DIGS;
				dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
				fp8_exp(c, a, d);
				TEST_ASSERT(fp8_cmp(c, b) == RLC_EQ, end);
			} TEST_END;

			TEST_BEGIN("frobenius and squared frobenius are consistent") {
				fp8_rand(a);
				fp8_frb(b, a, 1);
				fp8_frb(b, b, 1);
				fp8_frb(c, a, 2);
				TEST_ASSERT(fp8_cmp(c, b) == RLC_EQ, end);
			} TEST_END;
		}
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp8_free(a);
	fp8_free(b);
	fp8_free(c);
	bn_free(d);
	return code;
}


static int memory12(void) {
	err_t e;
	int code = RLC_ERR;
	fp12_t a;

	fp12_null(a);

	TRY {
		TEST_BEGIN("memory can be allocated") {
			fp12_new(a);
			fp12_free(a);
		} TEST_END;
	} CATCH(e) {
		switch (e) {
			case ERR_NO_MEMORY:
				util_print("FATAL ERROR!\n");
				ERROR(end);
				break;
		}
	}
	(void)a;
	code = RLC_OK;
  end:
	return code;
}

static int util12(void) {
	int code = RLC_ERR;
	uint8_t bin[12 * RLC_FP_BYTES];
	fp12_t a, b, c;
	dig_t d;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);

		TEST_BEGIN("comparison is consistent") {
			fp12_rand(a);
			fp12_rand(b);
			if (fp12_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp12_cmp(b, a) == RLC_NE, end);
			}
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_rand(c);
			if (fp12_cmp(a, c) != RLC_EQ) {
				fp12_copy(c, a);
				TEST_ASSERT(fp12_cmp(c, a) == RLC_EQ, end);
			}
			if (fp12_cmp(b, c) != RLC_EQ) {
				fp12_copy(c, b);
				TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
			}
		}
		TEST_END;

		TEST_BEGIN("negation is consistent") {
			fp12_rand(a);
			fp12_neg(b, a);
			if (fp12_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp12_cmp(b, a) == RLC_NE, end);
			}
			fp12_neg(b, b);
			TEST_ASSERT(fp12_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and comparison are consistent") {
			do {
				fp12_rand(a);
			} while (fp12_is_zero(a));
			fp12_zero(c);
			TEST_ASSERT(fp12_cmp(a, c) == RLC_NE, end);
			TEST_ASSERT(fp12_cmp(c, a) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to random and comparison are consistent") {
			do {
				fp12_rand(a);
			} while (fp12_is_zero(a));
			fp12_zero(c);
			TEST_ASSERT(fp12_cmp(a, c) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and zero test are consistent") {
			fp12_zero(a);
			TEST_ASSERT(fp12_is_zero(a), end);
		}
		TEST_END;

		TEST_BEGIN("assignment to a constant and comparison are consistent") {
			rand_bytes((uint8_t *)&d, (RLC_DIG / 8));
			fp12_set_dig(a, d);
			TEST_ASSERT(fp12_cmp_dig(a, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("reading and writing a finite field element are consistent") {
			fp12_rand(a);
			fp12_write_bin(bin, sizeof(bin), a, 0);
			fp12_read_bin(b, bin, sizeof(bin));
			TEST_ASSERT(fp12_cmp(a, b) == RLC_EQ, end);
			fp12_conv_cyc(a, a);
			fp12_write_bin(bin, 8 * RLC_FP_BYTES, a, 1);
			fp12_read_bin(b, bin, 8 * RLC_FP_BYTES);
			TEST_ASSERT(fp12_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("getting the size of a finite field element is correct") {
			fp12_rand(a);
			TEST_ASSERT(fp12_size_bin(a, 0) == 12 * RLC_FP_BYTES, end);
			fp12_conv_cyc(a, a);
			TEST_ASSERT(fp12_size_bin(a, 1) == 8 * RLC_FP_BYTES, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	return code;
}

static int addition12(void) {
	int code = RLC_ERR;
	fp12_t a, b, c, d, e;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);
	fp12_null(d);
	fp12_null(e);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);
		fp12_new(d);
		fp12_new(e);

		TEST_BEGIN("addition is commutative") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_add(d, a, b);
			fp12_add(e, b, a);
			TEST_ASSERT(fp12_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition is associative") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_rand(c);
			fp12_add(d, a, b);
			fp12_add(d, d, c);
			fp12_add(e, b, c);
			fp12_add(e, a, e);
			TEST_ASSERT(fp12_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has identity") {
			fp12_rand(a);
			fp12_zero(d);
			fp12_add(e, a, d);
			TEST_ASSERT(fp12_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has inverse") {
			fp12_rand(a);
			fp12_neg(d, a);
			fp12_add(e, a, d);
			TEST_ASSERT(fp12_is_zero(e), end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	fp12_free(d);
	fp12_free(e);
	return code;
}

static int subtraction12(void) {
	int code = RLC_ERR;
	fp12_t a, b, c, d;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);
	fp12_null(d);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);
		fp12_new(d);

		TEST_BEGIN("subtraction is anti-commutative") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_sub(c, a, b);
			fp12_sub(d, b, a);
			fp12_neg(d, d);
			TEST_ASSERT(fp12_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has identity") {
			fp12_rand(a);
			fp12_zero(c);
			fp12_sub(d, a, c);
			TEST_ASSERT(fp12_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has inverse") {
			fp12_rand(a);
			fp12_sub(c, a, a);
			TEST_ASSERT(fp12_is_zero(c), end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	fp12_free(d);
	return code;
}

static int multiplication12(void) {
	int code = RLC_ERR;
	fp12_t a, b, c, d, e, f;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);
	fp12_null(d);
	fp12_null(e);
	fp12_null(f);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);
		fp12_new(d);
		fp12_new(e);
		fp12_new(f);

		TEST_BEGIN("multiplication is commutative") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_mul(d, a, b);
			fp12_mul(e, b, a);
			TEST_ASSERT(fp12_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is associative") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_rand(c);
			fp12_mul(d, a, b);
			fp12_mul(d, d, c);
			fp12_mul(e, b, c);
			fp12_mul(e, a, e);
			TEST_ASSERT(fp12_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is distributive") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_rand(c);
			fp12_add(d, a, b);
			fp12_mul(d, c, d);
			fp12_mul(e, c, a);
			fp12_mul(f, c, b);
			fp12_add(e, e, f);
			TEST_ASSERT(fp12_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has identity") {
			fp12_zero(d);
			fp_set_dig(d[0][0][0], 1);
			fp12_mul(e, a, d);
			TEST_ASSERT(fp12_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has zero property") {
			fp12_zero(d);
			fp12_mul(e, a, d);
			TEST_ASSERT(fp12_is_zero(e), end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic multiplication is correct") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_mul(c, a, b);
			fp12_mul_basic(d, a, b);
			TEST_ASSERT(fp12_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced multiplication is correct") {
			fp12_rand(a);
			fp12_rand(b);
			fp12_mul(c, a, b);
			fp12_mul_lazyr(d, a, b);
			TEST_ASSERT(fp12_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("sparse multiplication is correct") {
			fp12_rand(a);
			fp12_rand(b);
			/* Test common case among configurations. */
			fp2_zero(b[0][1]);
			fp_zero(b[0][0][1]);
			fp2_zero(b[0][2]);
			fp2_zero(b[1][0]);
			fp_zero(b[1][1][1]);
			fp2_zero(b[1][2]);
			fp12_mul(c, a, b);
			fp12_mul_dxs(d, a, b);
			TEST_ASSERT(fp12_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic sparse multiplication is correct") {
			fp12_rand(a);
			fp12_rand(b);
			fp2_zero(b[0][1]);
			fp2_zero(b[0][2]);
			fp2_zero(b[1][2]);
			fp12_mul_dxs(c, a, b);
			fp12_mul_dxs_basic(d, a, b);
			TEST_ASSERT(fp12_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced sparse multiplication is correct") {
			fp12_rand(a);
			fp12_rand(b);
			fp2_zero(b[0][1]);
			fp2_zero(b[0][2]);
			fp2_zero(b[1][2]);
			fp12_mul_dxs(c, a, b);
			fp12_mul_dxs_lazyr(d, a, b);
			TEST_ASSERT(fp12_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	} CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	fp12_free(d);
	fp12_free(e);
	fp12_free(f);
	return code;
}

static int squaring12(void) {
	int code = RLC_ERR;
	fp12_t a, b, c;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);

		TEST_BEGIN("squaring is correct") {
			fp12_rand(a);
			fp12_mul(b, a, a);
			fp12_sqr(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic squaring is correct") {
			fp12_rand(a);
			fp12_sqr(b, a);
			fp12_sqr_basic(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced squaring is correct") {
			fp12_rand(a);
			fp12_sqr(b, a);
			fp12_sqr_lazyr(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	return code;
}

static int cyclotomic12(void) {
	int code = RLC_ERR;
	fp12_t a, b, c, d[2], e[2];
	bn_t f;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);
	fp12_null(d[0]);
	fp12_null(d[1])
	fp12_null(e[0]);
	bn_null(f);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);
		fp12_new(d[0]);
		fp12_new(d[1]);
		fp12_new(e[0]);
		fp12_new(e[1]);
		bn_new(f);

		TEST_BEGIN("cyclotomic test is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			TEST_ASSERT(fp12_test_cyc(a) == 1, end);
		} TEST_END;

		TEST_BEGIN("compression in cyclotomic subgroup is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp12_back_cyc(c, a);
			TEST_ASSERT(fp12_cmp(a, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("simultaneous compression in cyclotomic subgroup is correct") {
			fp12_rand(d[0]);
			fp12_rand(d[1]);
			fp12_conv_cyc(d[0], d[0]);
			fp12_conv_cyc(d[1], d[1]);
			fp12_back_cyc_sim(e, d, 2);
			TEST_ASSERT(fp12_cmp(d[0], e[0]) == RLC_EQ &&
					fp12_cmp(d[1], e[1]) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("cyclotomic squaring is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp12_sqr(b, a);
			fp12_sqr_cyc(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_BEGIN("basic cyclotomic squaring is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp12_sqr_cyc(b, a);
			fp12_sqr_cyc_basic(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced cyclotomic squaring is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp12_sqr_cyc(b, a);
			fp12_sqr_cyc_lazyr(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("compressed squaring is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp2_zero(b[0][0]);
			fp2_zero(b[1][1]);
			fp2_zero(c[0][0]);
			fp2_zero(c[1][1]);
			fp12_sqr(b, a);
			fp12_sqr_pck(c, a);
			fp12_back_cyc(c, c);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_BEGIN("basic compressed squaring is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp2_zero(b[0][0]);
			fp2_zero(b[1][1]);
			fp2_zero(c[0][0]);
			fp2_zero(c[1][1]);
			fp12_sqr_pck(b, a);
			fp12_sqr_pck_basic(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced compressed squaring is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp2_zero(b[0][0]);
			fp2_zero(b[1][1]);
			fp2_zero(c[0][0]);
			fp2_zero(c[1][1]);
			fp12_sqr_pck(b, a);
			fp12_sqr_pck_lazyr(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

        TEST_BEGIN("cyclotomic exponentiation is correct") {
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			bn_zero(f);
			fp12_exp_cyc(c, a, f);
			TEST_ASSERT(fp12_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(f, 1);
			fp12_exp_cyc(c, a, f);
			TEST_ASSERT(fp12_cmp(c, a) == RLC_EQ, end);
			bn_rand(f, RLC_POS, RLC_FP_BITS);
			fp12_exp(b, a, f);
			fp12_exp_cyc(c, a, f);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
			bn_rand(f, RLC_POS, RLC_FP_BITS);
			fp12_exp_cyc(b, a, f);
			bn_neg(f, f);
			fp12_exp_cyc(c, a, f);
			fp12_inv_uni(c, c);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
        } TEST_END;

		TEST_BEGIN("sparse cyclotomic exponentiation is correct") {
			int g[3] = {0, 0, RLC_FP_BITS - 1};
			do {
				bn_rand(f, RLC_POS, RLC_DIG);
				g[1] = f->dp[0] % RLC_FP_BITS;
			} while (g[1] == 0 || g[1] == RLC_FP_BITS - 1);
			bn_set_2b(f, RLC_FP_BITS - 1);
			bn_set_bit(f, g[1], 1);
			bn_set_bit(f, 0, 1);
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp12_exp(b, a, f);
			fp12_exp_cyc_sps(c, a, g, 3);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
			/*g[0] = 0;
			fp12_exp_cyc_sps(c, a, g, 0);
			TEST_ASSERT(fp12_cmp_dig(c, 1) == RLC_EQ, end);*/
			g[0] = 0;
			fp12_exp_cyc_sps(c, a, g, 1);
			TEST_ASSERT(fp12_cmp(c, a) == RLC_EQ, end);
			g[0] = -1;
			fp12_exp_cyc_sps(b, a, g, 1);
			fp12_inv(b, b);
			fp12_sqr_cyc(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	fp12_free(d[0]);
	fp12_free(d[1]);
	fp12_free(e[0]);
	fp12_free(e[1]);
	bn_free(f);
	return code;
}

static int inversion12(void) {
	int code = RLC_ERR;
	fp12_t a, b, c;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);

		TEST_BEGIN("inversion is correct") {
			do {
				fp12_rand(a);
			} while (fp12_is_zero(a));
			fp12_inv(b, a);
			fp12_mul(c, a, b);
			fp12_zero(b);
			fp_set_dig(b[0][0][0], 1);
			TEST_ASSERT(fp12_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("inversion of a unitary element is correct") {
			do {
				fp12_rand(a);
			} while (fp12_is_zero(a));
			fp12_conv_uni(a, a);
			fp12_inv(b, a);
			fp12_inv_uni(c, a);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	return code;
}

static int exponentiation12(void) {
	int code = RLC_ERR;
	fp12_t a, b, c;
	bn_t d;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);
	bn_null(d);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);
		bn_new(d);

		TEST_BEGIN("exponentiation is correct") {
			fp12_rand(a);
			bn_zero(d);
			fp12_exp(c, a, d);
			TEST_ASSERT(fp12_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(d, 1);
			fp12_exp(c, a, d);
			TEST_ASSERT(fp12_cmp(c, a) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_FP_BITS);
			fp12_exp(b, a, d);
			bn_neg(d, d);
			fp12_exp(c, a, d);
			fp12_inv(c, c);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_DIG);
			fp12_exp_dig(b, a, d->dp[0]);
			fp12_exp(c, a, d);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and exponentiation are consistent") {
			fp12_rand(a);
			fp12_frb(b, a, 0);
			TEST_ASSERT(fp12_cmp(a, b) == RLC_EQ, end);
			fp12_frb(b, a, 1);
			d->sign = RLC_POS;
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			fp12_exp(c, a, d);
			TEST_ASSERT(fp12_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and squared frobenius are consistent") {
			fp12_rand(a);
			fp12_frb(b, a, 1);
			fp12_frb(b, b, 1);
			fp12_frb(c, a, 2);
			TEST_ASSERT(fp12_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and cubed frobenius are consistent") {
			fp12_rand(a);
			fp12_frb(b, a, 1);
			fp12_frb(b, b, 1);
			fp12_frb(b, b, 1);
			fp12_frb(c, a, 3);
			TEST_ASSERT(fp12_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	bn_free(d);
	return code;
}

static int compression12(void) {
	int code = RLC_ERR;
	fp12_t a, b, c;

	fp12_null(a);
	fp12_null(b);
	fp12_null(c);

	TRY {
		fp12_new(a);
		fp12_new(b);
		fp12_new(c);

		TEST_BEGIN("compression is consistent") {
			fp12_rand(a);
			fp12_pck(b, a);
			TEST_ASSERT(fp12_upk(c, b) == 1, end);
			TEST_ASSERT(fp12_cmp(a, c) == RLC_EQ, end);
			fp12_rand(a);
			fp12_conv_cyc(a, a);
			fp12_pck(b, a);
			TEST_ASSERT(fp12_upk(c, b) == 1, end);
			TEST_ASSERT(fp12_cmp(a, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp12_free(a);
	fp12_free(b);
	fp12_free(c);
	return code;
}

static int memory18(void) {
	err_t e;
	int code = RLC_ERR;
	fp18_t a;

	fp18_null(a);

	TRY {
		TEST_BEGIN("memory can be allocated") {
			fp18_new(a);
			fp18_free(a);
		} TEST_END;
	} CATCH(e) {
		switch (e) {
			case ERR_NO_MEMORY:
				util_print("FATAL ERROR!\n");
				ERROR(end);
				break;
		}
	}
	(void)a;
	code = RLC_OK;
  end:
	return code;
}

static int util18(void) {
	int code = RLC_ERR;
	fp18_t a, b, c;
	dig_t d;

	fp18_null(a);
	fp18_null(b);
	fp18_null(c);

	TRY {
		fp18_new(a);
		fp18_new(b);
		fp18_new(c);

		TEST_BEGIN("comparison is consistent") {
			fp18_rand(a);
			fp18_rand(b);
			if (fp18_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp18_cmp(b, a) == RLC_NE, end);
			}
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_rand(c);
			if (fp18_cmp(a, c) != RLC_EQ) {
				fp18_copy(c, a);
				TEST_ASSERT(fp18_cmp(c, a) == RLC_EQ, end);
			}
			if (fp18_cmp(b, c) != RLC_EQ) {
				fp18_copy(c, b);
				TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
			}
		}
		TEST_END;

		TEST_BEGIN("negation is consistent") {
			fp18_rand(a);
			fp18_neg(b, a);
			if (fp18_cmp(a, b) != RLC_EQ) {
				TEST_ASSERT(fp18_cmp(b, a) == RLC_NE, end);
			}
			fp18_neg(b, b);
			TEST_ASSERT(fp18_cmp(a, b) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and comparison are consistent") {
			do {
				fp18_rand(a);
			} while (fp18_is_zero(a));
			fp18_zero(c);
			TEST_ASSERT(fp18_cmp(a, c) == RLC_NE, end);
			TEST_ASSERT(fp18_cmp(c, a) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and comparison are consistent") {
			do {
				fp18_rand(a);
			} while (fp18_is_zero(a));
			fp18_zero(c);
			TEST_ASSERT(fp18_cmp(a, c) == RLC_NE, end);
			TEST_ASSERT(fp18_cmp(c, a) == RLC_NE, end);
		}
		TEST_END;


		TEST_BEGIN("assignment to random and comparison are consistent") {
			do {
				fp18_rand(a);
			} while (fp18_is_zero(a));
			fp18_zero(c);
			TEST_ASSERT(fp18_cmp(a, c) == RLC_NE, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to zero and zero test are consistent") {
			fp18_zero(a);
			TEST_ASSERT(fp18_is_zero(a), end);
		}
		TEST_END;

		TEST_BEGIN("assignment to a constant and comparison are consistent") {
			rand_bytes((uint8_t *)&d, (RLC_DIG / 8));
			fp18_set_dig(a, d);
			TEST_ASSERT(fp18_cmp_dig(a, d) == RLC_EQ, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp18_free(a);
	fp18_free(b);
	fp18_free(c);
	return code;
}

static int addition18(void) {
	int code = RLC_ERR;
	fp18_t a, b, c, d, e;

	fp18_null(a);
	fp18_null(b);
	fp18_null(c);
	fp18_null(d);
	fp18_null(e);

	TRY {
		fp18_new(a);
		fp18_new(b);
		fp18_new(c);
		fp18_new(d);
		fp18_new(e);

		TEST_BEGIN("addition is commutative") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_add(d, a, b);
			fp18_add(e, b, a);
			TEST_ASSERT(fp18_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition is associative") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_rand(c);
			fp18_add(d, a, b);
			fp18_add(d, d, c);
			fp18_add(e, b, c);
			fp18_add(e, a, e);
			TEST_ASSERT(fp18_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has identity") {
			fp18_rand(a);
			fp18_zero(d);
			fp18_add(e, a, d);
			TEST_ASSERT(fp18_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("addition has inverse") {
			fp18_rand(a);
			fp18_neg(d, a);
			fp18_add(e, a, d);
			TEST_ASSERT(fp18_is_zero(e), end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp18_free(a);
	fp18_free(b);
	fp18_free(c);
	fp18_free(d);
	fp18_free(e);
	return code;
}

static int subtraction18(void) {
	int code = RLC_ERR;
	fp18_t a, b, c, d;

	TRY {
		fp18_new(a);
		fp18_new(b);
		fp18_new(c);
		fp18_new(d);

		TEST_BEGIN("subtraction is anti-commutative") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_sub(c, a, b);
			fp18_sub(d, b, a);
			fp18_neg(d, d);
			TEST_ASSERT(fp18_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has identity") {
			fp18_rand(a);
			fp18_zero(c);
			fp18_sub(d, a, c);
			TEST_ASSERT(fp18_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("subtraction has inverse") {
			fp18_rand(a);
			fp18_sub(c, a, a);
			TEST_ASSERT(fp18_is_zero(c), end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp18_free(a);
	fp18_free(b);
	fp18_free(c);
	fp18_free(d);
	return code;
}

static int multiplication18(void) {
	int code = RLC_ERR;
	fp18_t a, b, c, d, e, f;

	TRY {
		fp18_new(a);
		fp18_new(b);
		fp18_new(c);
		fp18_new(d);
		fp18_new(e);
		fp18_new(f);

		TEST_BEGIN("multiplication is commutative") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_mul(d, a, b);
			fp18_mul(e, b, a);
			TEST_ASSERT(fp18_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is associative") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_rand(c);
			fp18_mul(d, a, b);
			fp18_mul(d, d, c);
			fp18_mul(e, b, c);
			fp18_mul(e, a, e);
			TEST_ASSERT(fp18_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication is distributive") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_rand(c);
			fp18_add(d, a, b);
			fp18_mul(d, c, d);
			fp18_mul(e, c, a);
			fp18_mul(f, c, b);
			fp18_add(e, e, f);
			TEST_ASSERT(fp18_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has identity") {
			fp18_zero(d);
			fp_set_dig(d[0][0][0], 1);
			fp18_mul(e, a, d);
			TEST_ASSERT(fp18_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multiplication has zero property") {
			fp18_zero(d);
			fp18_mul(e, a, d);
			TEST_ASSERT(fp18_is_zero(e), end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic multiplication is correct") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_mul(c, a, b);
			fp18_mul_basic(d, a, b);
			TEST_ASSERT(fp18_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced multiplication is correct") {
			fp18_rand(a);
			fp18_rand(b);
			fp18_mul(c, a, b);
			fp18_mul_lazyr(d, a, b);
			TEST_ASSERT(fp18_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	} CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp18_free(a);
	fp18_free(b);
	fp18_free(c);
	fp18_free(d);
	fp18_free(e);
	fp18_free(f);
	return code;
}

static int squaring18(void) {
	int code = RLC_ERR;
	fp18_t a, b, c;

	fp18_null(a);
	fp18_null(b);
	fp18_null(c);

	TRY {
		fp18_new(a);
		fp18_new(b);
		fp18_new(c);

		TEST_BEGIN("squaring is correct") {
			fp18_rand(a);
			fp18_mul(b, a, a);
			fp18_sqr(c, a);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC | !defined(STRIP)
		TEST_BEGIN("basic squaring is correct") {
			fp18_rand(a);
			fp18_sqr(b, a);
			fp18_sqr_basic(c, a);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy reduced squaring is correct") {
			fp18_rand(a);
			fp18_sqr(b, a);
			fp18_sqr_lazyr(c, a);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp18_free(a);
	fp18_free(b);
	fp18_free(c);
	return code;
}

static int cyclotomic18(void) {
	int code = RLC_ERR;
	fp18_t a, b, c, d[2], e[2];
	bn_t f;

	fp18_null(a);
	fp18_null(b);
	fp18_null(c);
	fp18_null(d[0]);
	fp18_null(d[1]);
	fp18_null(e[0]);
	fp18_null(e[1]);
	bn_null(f);

	TRY {
		fp18_new(a);
		fp18_new(b);
		fp18_new(c);
		fp18_new(d[0]);
		fp18_new(d[1]);
		fp18_new(e[0]);
		fp18_new(e[1]);
		bn_new(f);

		TEST_BEGIN("cyclotomic test is correct") {
			fp18_rand(a);
			fp18_conv_cyc(a, a);
			TEST_ASSERT(fp18_test_cyc(a) == 1, end);
		} TEST_END;

		TEST_BEGIN("compression in cyclotomic subgroup is correct") {
			fp18_rand(a);
			fp18_conv_cyc(a, a);
			fp18_back_cyc(c, a);
			TEST_ASSERT(fp18_cmp(a, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("simultaneous compression in cyclotomic subgroup is correct") {
			fp18_rand(d[0]);
			fp18_rand(d[1]);
			fp18_conv_cyc(d[0], d[0]);
			fp18_conv_cyc(d[1], d[1]);
			fp18_back_cyc_sim(e, d, 2);
			TEST_ASSERT(fp18_cmp(d[0], e[0]) == RLC_EQ &&
					fp18_cmp(d[1], e[1]) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("cyclotomic squaring is correct") {
			fp18_rand(a);
			fp18_conv_cyc(a, a);
			fp18_sqr(b, a);
			fp18_sqr_cyc(c, a);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("compressed squaring is correct") {
			fp18_rand(a);
			fp18_conv_cyc(a, a);
			fp2_zero(b[0][0]);
			fp2_zero(b[1][1]);
			fp2_zero(c[0][0]);
			fp2_zero(c[1][1]);
			fp18_sqr(b, a);
			fp18_sqr_pck(c, a);
			fp18_back_cyc(c, c);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

        TEST_BEGIN("cyclotomic exponentiation is correct") {
			bn_rand(f, RLC_POS, RLC_FP_BITS);
			fp18_rand(a);
			fp18_conv_cyc(a, a);
			fp18_exp(b, a, f);
			fp18_exp_cyc(c, a, f);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
			bn_zero(f);
			fp18_exp_cyc(c, a, f);
			TEST_ASSERT(fp18_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(f, 1);
			fp18_exp_cyc(c, a, f);
			TEST_ASSERT(fp18_cmp(c, a) == RLC_EQ, end);
			bn_rand(f, RLC_POS, RLC_FP_BITS);
			fp18_exp_cyc(b, a, f);
			bn_neg(f, f);
			fp18_exp_cyc(c, a, f);
			fp18_inv(c, c);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
        } TEST_END;

		TEST_BEGIN("sparse cyclotomic exponentiation is correct") {
			int g[3] = {0, 0, RLC_FP_BITS - 1};
			do {
				bn_rand(f, RLC_POS, RLC_DIG);
				g[1] = f->dp[0] % RLC_FP_BITS;
			} while (g[1] == 0 || g[1] == RLC_FP_BITS - 1);
			bn_set_2b(f, RLC_FP_BITS - 1);
			bn_set_bit(f, g[1], 1);
			bn_set_bit(f, 0, 1);
			fp18_rand(a);
			fp18_conv_cyc(a, a);
			fp18_exp(b, a, f);
			fp18_exp_cyc_sps(c, a, g, 3);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
			g[0] = 0;
			fp18_exp_cyc_sps(c, a, g, 0);
			TEST_ASSERT(fp18_cmp_dig(c, 1) == RLC_EQ, end);
			g[0] = 0;
			fp18_exp_cyc_sps(c, a, g, 1);
			TEST_ASSERT(fp18_cmp(c, a) == RLC_EQ, end);
			g[0] = -1;
			fp18_exp_cyc_sps(b, a, g, 1);
			fp18_inv(b, b);
			fp18_sqr_cyc(c, a);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp18_free(a);
	fp18_free(b);
	fp18_free(c);
	fp18_free(d[0]);
	fp18_free(d[1]);
	fp18_free(e[0]);
	fp18_free(e[1]);
	bn_free(f);
	return code;
}

static int inversion18(void) {
	int code = RLC_ERR;
	fp18_t a, b, c;

	fp18_null(a);
	fp18_null(b);
	fp18_null(c);

	TRY {
		fp18_new(a);
		fp18_new(b);
		fp18_new(c);

		TEST_BEGIN("inversion is correct") {
			do {
				fp18_rand(a);
			} while (fp18_is_zero(a));
			fp18_inv(b, a);
			fp18_mul(c, a, b);
			fp18_zero(b);
			fp_set_dig(b[0][0][0], 1);
			TEST_ASSERT(fp18_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("inversion of a unitary element is correct") {
			do {
				fp18_rand(a);
			} while (fp18_is_zero(a));
			fp18_conv_uni(a, a);
			fp18_inv(b, a);
			fp18_inv_uni(c, a);
			TEST_ASSERT(fp18_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp18_free(a);
	fp18_free(b);
	fp18_free(c);
	return code;
}

static int exponentiation18(void) {
	int code = RLC_ERR;
	fp18_t a, b, c;
	bn_t d;

	fp18_null(a);
	fp18_null(b);
	fp18_null(c);
	bn_null(d);

	TRY {
		fp18_new(a);
		fp18_new(b);
		fp18_new(c);
		bn_new(d);

		TEST_BEGIN("exponentiation is correct") {
			fp18_rand(a);
			bn_zero(d);
			fp18_exp(c, a, d);
			TEST_ASSERT(fp18_cmp_dig(c, 1) == RLC_EQ, end);
			bn_set_dig(d, 1);
			fp18_exp(c, a, d);
			TEST_ASSERT(fp18_cmp(c, a) == RLC_EQ, end);
			bn_rand(d, RLC_POS, RLC_FP_BITS);
			fp18_exp(b, a, d);
			bn_neg(d, d);
			fp18_exp(c, a, d);
			fp18_inv(c, c);
			TEST_ASSERT(fp12_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and exponentiation are consistent") {
			fp18_rand(a);
			fp18_frb(b, a, 0);
			TEST_ASSERT(fp18_cmp(a, b) == RLC_EQ, end);
			fp18_frb(b, a, 1);
			d->sign = RLC_POS;
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			fp18_exp(c, a, d);
			TEST_ASSERT(fp18_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and squared frobenius are consistent") {
			fp18_rand(a);
			fp18_frb(b, a, 1);
			fp18_frb(b, b, 1);
			fp18_frb(c, a, 2);
			TEST_ASSERT(fp18_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and cubed frobenius are consistent") {
			fp18_rand(a);
			fp18_frb(b, a, 1);
			fp18_frb(b, b, 1);
			fp18_frb(b, b, 1);
			fp18_frb(c, a, 3);
			TEST_ASSERT(fp18_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and quartered frobenius are consistent") {
			fp18_rand(a);
			fp18_frb(b, a, 1);
			fp18_frb(b, b, 1);
			fp18_frb(b, b, 1);
			fp18_frb(b, b, 1);
			fp18_frb(c, a, 4);
			TEST_ASSERT(fp18_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("frobenius and quintered frobenius are consistent") {
			fp18_rand(a);
			fp18_frb(b, a, 1);
			fp18_frb(b, b, 1);
			fp18_frb(b, b, 1);
			fp18_frb(b, b, 1);
			fp18_frb(b, b, 1);
			fp18_frb(c, a, 5);
			TEST_ASSERT(fp18_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	fp18_free(a);
	fp18_free(b);
	fp18_free(c);
	bn_free(d);
	return code;
}

int main(void) {
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	util_banner("Tests for the FPX module", 0);

	/* Try using a pairing-friendly curve for faster exponentiation method. */
	if (pc_param_set_any() != RLC_OK) {
		//If it does not work, try a tower-friendly field. 
		if (fp_param_set_any_tower() == RLC_ERR) {
			THROW(ERR_NO_FIELD);
			core_clean();
			return 0;
		}
	}

	/* Only execute these if there is an assigned quadratic non-residue. */
	if (fp_prime_get_qnr()) {
		util_print("\n-- Quadratic extension: %d as QNR\n", fp_prime_get_qnr());
		util_banner("Utilities:", 1);

		if (memory2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (util2() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Arithmetic:", 1);

		if (addition2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (subtraction2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (multiplication2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (squaring2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (inversion2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (exponentiation2() != RLC_OK) {
			core_clean();
			return 1;
		}

#ifdef FP_QNRES
		/* Restrict compression to p = 3 mod 4 for the moment. */
		if (compression2() != RLC_OK) {
			core_clean();
			return 1;
		}
#endif

		if (square_root2() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	/* Only execute these if there is an assigned cubic non-residue. */
	if (fp_prime_get_cnr()) {
		util_print("\n-- Cubic extension: %d as CNR\n", fp_prime_get_cnr());
		util_banner("Utilities:", 1);

		if (memory3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (util3() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Arithmetic:", 1);

		if (addition3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (subtraction3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (multiplication3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (squaring3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (inversion3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (exponentiation3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (square_root3() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	/* Fp^4 is defined as a quadratic extension of Fp^2. */
	if (fp_prime_get_qnr()) {
		util_print("\n-- Quartic extension: (i + %d) as QNR\n",
				fp2_field_get_qnr());
		util_banner("Utilities:", 1);

		if (memory4() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (util4() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Arithmetic:", 1);

		if (addition4() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (subtraction4() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling4() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (multiplication4() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (squaring4() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (inversion4() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (exponentiation4() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	/* Fp^6 is defined as a cubic extension of Fp^2. */
	if (fp_prime_get_qnr()) {
		util_print("\n-- Sextic extension: (i + %d) as QNR\n",
				fp2_field_get_qnr());
		util_banner("Utilities:", 1);

		if (memory6() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (util6() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Arithmetic:", 1);

		if (addition6() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (subtraction6() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling6() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (multiplication6() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (squaring6() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (inversion6() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (exponentiation6() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Octic extension:", 0);
		util_banner("Utilities:", 1);

		if (memory8() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (util8() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Arithmetic:", 1);

		if (addition8() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (subtraction8() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling8() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (multiplication8() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (squaring8() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (cyclotomic8() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (inversion8() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (exponentiation8() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Dodecic extension:", 0);
		util_banner("Utilities:", 1);

		if (memory12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (util12() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Arithmetic:", 1);

		if (addition12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (subtraction12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (multiplication12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (squaring12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (cyclotomic12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (inversion12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (exponentiation12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (compression12() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	/*
	 * Fp^18 is defined as a cubic extension of Fp^6, which in turn is defined
	 * as a quadratic extension of Fp^3, but defined as a cubic extension of
	 * Fp^2. Only run these if there is a cubic and quadratic non-residue in Fp.
	 */
	if (fp_prime_get_cnr() == fp_prime_get_qnr()) {

		util_banner("Octdecic extension:", 0);
		util_banner("Utilities:", 1);

		if (memory18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (util18() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Arithmetic:", 1);

		if (addition18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (subtraction18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (multiplication18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (squaring18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (cyclotomic18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (inversion18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (exponentiation18() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	util_banner("All tests have passed.\n", 0);

	core_clean();
	return 0;
}
