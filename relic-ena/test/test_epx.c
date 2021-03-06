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
 * Tests for elliptic curves defined over extensions of prime fields.
 *
 * @ingroup test
 */

#include <stdio.h>

#include "relic.h"
#include "relic_test.h"

#if FP_PRIME==508
static int memory(void) {
	err_t e;
	int code = RLC_ERR;
	ep3_t a;

	ep3_null(a);

	TRY {
		TEST_BEGIN("memory can be allocated") {
			ep3_new(a);
			ep3_free(a);
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
#else
static int memory(void) {
	err_t e;
	int code = RLC_ERR;
	ep2_t a;

	ep2_null(a);

	TRY{
		TEST_BEGIN("memory can be allocated") {
			ep2_new(a);
			ep2_free(a);
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
#endif

int util(void) {
	int l, code = RLC_ERR;
	ep2_t a, b, c;
	bn_t n;
	uint8_t bin[4 * RLC_FP_BYTES + 1];

	ep2_null(a);
	ep2_null(b);
	ep2_null(c);
	bn_null(n);

	TRY {
		ep2_new(a);
		ep2_new(b);
		ep2_new(c);
		bn_new(n);

		TEST_BEGIN("comparison is consistent") {
			ep2_rand(a);
			ep2_rand(b);
			TEST_ASSERT(ep2_cmp(a, b) != RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_rand(c);
			/* Compare points in affine coordinates. */
			if (ep2_cmp(a, c) != RLC_EQ) {
				ep2_copy(c, a);
				TEST_ASSERT(ep2_cmp(c, a) == RLC_EQ, end);
			}
			if (ep2_cmp(b, c) != RLC_EQ) {
				ep2_copy(c, b);
				TEST_ASSERT(ep2_cmp(b, c) == RLC_EQ, end);
			}
			/* Compare with one point in projective. */
			ep2_dbl(c, a);
			ep2_norm(c, c);
			ep2_dbl(a, a);
			TEST_ASSERT(ep2_cmp(c, a) == RLC_EQ, end);
			TEST_ASSERT(ep2_cmp(a, c) == RLC_EQ, end);
			/* Compare with two points in projective. */
			ep2_dbl(c, c);
			ep2_dbl(a, a);
			TEST_ASSERT(ep2_cmp(c, a) == RLC_EQ, end);
			TEST_ASSERT(ep2_cmp(a, c) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("negation and comparison are consistent") {
			ep2_rand(a);
			ep2_neg(b, a);
			TEST_ASSERT(ep2_cmp(a, b) != RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN
				("assignment to random/infinity and comparison are consistent")
		{
			ep2_rand(a);
			ep2_set_infty(c);
			TEST_ASSERT(ep2_cmp(a, c) != RLC_EQ, end);
			TEST_ASSERT(ep2_cmp(c, a) != RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to infinity and infinity test are consistent") {
			ep2_set_infty(a);
			TEST_ASSERT(ep2_is_infty(a), end);
		}
		TEST_END;

		TEST_BEGIN("validity test is correct") {
			ep2_set_infty(a);
			TEST_ASSERT(ep2_is_valid(a), end);
			ep2_rand(a);
			TEST_ASSERT(ep2_is_valid(a), end);
			fp2_rand(a->x);
			TEST_ASSERT(!ep2_is_valid(a), end);
		}
		TEST_END;

		TEST_BEGIN("reading and writing a point are consistent") {
			for (int j = 0; j < 2; j++) {
				ep2_set_infty(a);
				l = ep2_size_bin(a, j);
				ep2_write_bin(bin, l, a, j);
				ep2_read_bin(b, bin, l);
				TEST_ASSERT(ep2_cmp(a, b) == RLC_EQ, end);
				ep2_rand(a);
				l = ep2_size_bin(a, j);
				ep2_write_bin(bin, l, a, j);
				ep2_read_bin(b, bin, l);
				TEST_ASSERT(ep2_cmp(a, b) == RLC_EQ, end);
				ep2_rand(a);
				ep2_dbl(a, a);
				l = ep2_size_bin(a, j);
				ep2_norm(a, a);
				ep2_write_bin(bin, l, a, j);
				ep2_read_bin(b, bin, l);
				TEST_ASSERT(ep2_cmp(a, b) == RLC_EQ, end);
			}
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	ep2_free(a);
	ep2_free(b);
	ep2_free(c);
	bn_free(n);
	return code;
}

int addition(void) {
	int code = RLC_ERR;
	ep2_t a, b, c, d, e;

	ep2_null(a);
	ep2_null(b);
	ep2_null(c);
	ep2_null(d);
	ep2_null(e);

	TRY {
		ep2_new(a);
		ep2_new(b);
		ep2_new(c);
		ep2_new(d);
		ep2_new(e);

		TEST_BEGIN("point addition is commutative") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_add(d, a, b);
			ep2_add(e, b, a);
			TEST_ASSERT(ep2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point addition is associative") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_rand(c);
			ep2_add(d, a, b);
			ep2_add(d, d, c);
			ep2_add(e, b, c);
			ep2_add(e, e, a);
			TEST_ASSERT(ep2_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point addition has identity") {
			ep2_rand(a);
			ep2_set_infty(d);
			ep2_add(e, a, d);
			TEST_ASSERT(ep2_cmp(e, a) == RLC_EQ, end);
			ep2_add(e, d, a);
			TEST_ASSERT(ep2_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point addition has inverse") {
			ep2_rand(a);
			ep2_neg(d, a);
			ep2_add(e, a, d);
			TEST_ASSERT(ep2_is_infty(e), end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("point addition in affine coordinates is correct") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_add(d, a, b);
			ep2_norm(d, d);
			ep2_add_basic(e, a, b);
			TEST_ASSERT(ep2_cmp(e, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
#if !defined(EP_MIXED) || !defined(STRIP)
		TEST_BEGIN("point addition in projective coordinates is correct") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_add_projc(a, a, b);
			ep2_rand(b);
			ep2_rand(c);
			ep2_add_projc(b, b, c);
			ep2_norm(b, b);
			/* a and b in projective coordinates. */
			ep2_add_projc(d, a, b);
			ep2_norm(d, d);
			ep2_norm(a, a);
			ep2_norm(b, b);
			ep2_add(e, a, b);
			ep2_norm(e, e);
			TEST_ASSERT(ep2_cmp(e, d) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("point addition in mixed coordinates (z2 = 1) is correct") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_add_projc(a, a, b);
			ep2_rand(b);
			/* a and b in projective coordinates. */
			ep2_add_projc(d, a, b);
			ep2_norm(d, d);
			/* a in affine coordinates. */
			ep2_norm(a, a);
			ep2_add(e, a, b);
			ep2_norm(e, e);
			TEST_ASSERT(ep2_cmp(e, d) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point addition in mixed coordinates (z1,z2 = 1) is correct") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_norm(a, a);
			ep2_norm(b, b);
			/* a and b in affine coordinates. */
			ep2_add(d, a, b);
			ep2_norm(d, d);
			ep2_add_projc(e, a, b);
			ep2_norm(e, e);
			TEST_ASSERT(ep2_cmp(e, d) == RLC_EQ, end);
		} TEST_END;
#endif

	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	ep2_free(a);
	ep2_free(b);
	ep2_free(c);
	ep2_free(d);
	ep2_free(e);
	return code;
}

int subtraction(void) {
	int code = RLC_ERR;
	ep2_t a, b, c, d;

	ep2_null(a);
	ep2_null(b);
	ep2_null(c);
	ep2_null(d);

	TRY {
		ep2_new(a);
		ep2_new(b);
		ep2_new(c);
		ep2_new(d);

		TEST_BEGIN("point subtraction is anti-commutative") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_sub(c, a, b);
			ep2_sub(d, b, a);
			ep2_neg(d, d);
			TEST_ASSERT(ep2_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("point subtraction has identity") {
			ep2_rand(a);
			ep2_set_infty(c);
			ep2_sub(d, a, c);
			TEST_ASSERT(ep2_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("point subtraction has inverse") {
			ep2_rand(a);
			ep2_sub(c, a, a);
			TEST_ASSERT(ep2_is_infty(c), end);
		}
		TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("point subtraction in affine coordinates is correct") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_sub(c, a, b);
			ep2_norm(c, c);
			ep2_sub_basic(d, a, b);
			TEST_ASSERT(ep2_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
#if !defined(EP_MIXED) || !defined(STRIP)
		TEST_BEGIN("point subtraction in projective coordinates is correct") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_add_projc(a, a, b);
			ep2_rand(b);
			ep2_rand(c);
			ep2_add_projc(b, b, c);
			/* a and b in projective coordinates. */
			ep2_sub_projc(c, a, b);
			ep2_norm(c, c);
			ep2_norm(a, a);
			ep2_norm(b, b);
			ep2_sub(d, a, b);
			ep2_norm(d, d);
			TEST_ASSERT(ep2_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("point subtraction in mixed coordinates (z2 = 1) is correct") {
			ep2_rand(a);
			ep2_rand(b);
			ep2_add_projc(a, a, b);
			ep2_rand(b);
			/* a and b in projective coordinates. */
			ep2_sub_projc(c, a, b);
			ep2_norm(c, c);
			/* a in affine coordinates. */
			ep2_norm(a, a);
			ep2_sub(d, a, b);
			ep2_norm(d, d);
			TEST_ASSERT(ep2_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN
				("point subtraction in mixed coordinates (z1,z2 = 1) is correct")
		{
			ep2_rand(a);
			ep2_rand(b);
			ep2_norm(a, a);
			ep2_norm(b, b);
			/* a and b in affine coordinates. */
			ep2_sub(c, a, b);
			ep2_norm(c, c);
			ep2_sub_projc(d, a, b);
			ep2_norm(d, d);
			TEST_ASSERT(ep2_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	ep2_free(a);
	ep2_free(b);
	ep2_free(c);
	ep2_free(d);
	return code;
}

int doubling(void) {
	int code = RLC_ERR;
	ep2_t a, b, c;

	ep2_null(a);
	ep2_null(b);
	ep2_null(c);

	TRY {
		ep2_new(a);
		ep2_new(b);
		ep2_new(c);

		TEST_BEGIN("point doubling is correct") {
			ep2_rand(a);
			ep2_add(b, a, a);
			ep2_dbl(c, a);
			TEST_ASSERT(ep2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("point doubling in affine coordinates is correct") {
			ep2_rand(a);
			ep2_dbl(b, a);
			ep2_norm(b, b);
			ep2_dbl_basic(c, a);
			TEST_ASSERT(ep2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
		TEST_BEGIN("point doubling in projective coordinates is correct") {
			ep2_rand(a);
			ep2_dbl_projc(a, a);
			/* a in projective coordinates. */
			ep2_dbl_projc(b, a);
			ep2_norm(b, b);
			ep2_norm(a, a);
			ep2_dbl(c, a);
			ep2_norm(c, c);
			TEST_ASSERT(ep2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point doubling in mixed coordinates (z1 = 1) is correct") {
			ep2_rand(a);
			ep2_dbl_projc(b, a);
			ep2_norm(b, b);
			ep2_dbl(c, a);
			ep2_norm(c, c);
			TEST_ASSERT(ep2_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	ep2_free(a);
	ep2_free(b);
	ep2_free(c);
	return code;
}

static int multiplication(void) {
	int code = RLC_ERR;
	bn_t n, k;
	ep2_t p, q, r;

	bn_null(n);
	bn_null(k);
	ep2_null(p);
	ep2_null(q);
	ep2_null(r);

	TRY {
		bn_new(n);
		bn_new(k);
		ep2_new(p);
		ep2_new(q);
		ep2_new(r);

		ep2_curve_get_gen(p);
		ep2_curve_get_ord(n);

		TEST_BEGIN("generator has the right order") {
			TEST_ASSERT(ep2_is_valid(p), end);
			ep2_mul(r, p, n);
			TEST_ASSERT(ep2_is_infty(r) == 1, end);
		} TEST_END;

		TEST_BEGIN("generator multiplication is correct") {
			bn_zero(k);
			ep2_mul_gen(r, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_gen(r, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep2_mul(q, p, k);
			ep2_mul_gen(r, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_gen(r, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;

#if EP_MUL == BASIC || !defined(STRIP)
		TEST_BEGIN("binary point multiplication is correct") {
			bn_zero(k);
			ep2_mul_basic(r, p, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_basic(r, p, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			ep2_rand(p);
			ep2_mul(r, p, n);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_rand_mod(k, n);
			ep2_mul(q, p, k);
			ep2_mul_basic(r, p, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_basic(r, p, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_MUL == MONTY || !defined(STRIP)
		TEST_BEGIN("sliding window point multiplication is correct") {
			bn_zero(k);
			ep2_mul_slide(r, p, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_slide(r, p, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			ep2_rand(p);
			ep2_mul(r, p, n);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_rand_mod(k, n);
			ep2_mul(q, p, k);
			ep2_mul_slide(r, p, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_slide(r, p, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		}
		TEST_END;
#endif

#if EP_MUL == MONTY || !defined(STRIP)
		TEST_BEGIN("montgomery laddering point multiplication is correct") {
			bn_zero(k);
			ep2_mul_monty(r, p, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_monty(r, p, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			ep2_rand(p);
			ep2_mul(r, p, n);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_rand_mod(k, n);
			ep2_mul(q, p, k);
			ep2_mul_monty(r, p, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_monty(r, p, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		}
		TEST_END;
#endif

#if EP_MUL == LWNAF || !defined(STRIP)
		TEST_BEGIN("left-to-right w-naf point multiplication is correct") {
			bn_zero(k);
			ep2_mul_lwnaf(r, p, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_lwnaf(r, p, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			ep2_rand(p);
			ep2_mul(r, p, n);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_rand_mod(k, n);
			ep2_mul(q, p, k);
			ep2_mul_lwnaf(r, p, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_lwnaf(r, p, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		}
		TEST_END;
#endif

		TEST_BEGIN("multiplication by digit is correct") {
			ep2_mul_dig(r, p, 0);
			TEST_ASSERT(ep2_is_infty(r), end);
			ep2_mul_dig(r, p, 1);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			bn_rand(k, RLC_POS, RLC_DIG);
			ep2_mul(q, p, k);
			ep2_mul_dig(r, p, k->dp[0]);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep2_free(p);
	ep2_free(q);
	ep2_free(r);
	return code;
}

static int fixed(void) {
	int code = RLC_ERR;
	bn_t n, k;
	ep2_t p, q, r, t[RLC_EPX_TABLE_MAX];

	bn_null(n);
	bn_null(k);
	ep2_null(p);
	ep2_null(q);
	ep2_null(r);

	for (int i = 0; i < RLC_EPX_TABLE_MAX; i++) {
		ep2_null(t[i]);
	}

	TRY {
		bn_new(n);
		bn_new(k);
		ep2_new(p);
		ep2_new(q);
		ep2_new(r);

		ep2_curve_get_gen(p);
		ep2_curve_get_ord(n);

		for (int i = 0; i < RLC_EP_TABLE; i++) {
			ep2_new(t[i]);
		}
		TEST_BEGIN("fixed point multiplication is correct") {
			ep2_rand(p);
			ep2_mul_pre(t, p);
			bn_zero(k);
			ep2_mul_fix(r, t, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_fix(r, t, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep2_mul(q, p, k);
			ep2_mul_fix(q, t, k);
			ep2_mul(r, p, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_fix(r, t, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE; i++) {
			ep2_free(t[i]);
		}

#if EP_FIX == BASIC || !defined(STRIP)
		for (int i = 0; i < RLC_EP_TABLE_BASIC; i++) {
			ep2_new(t[i]);
		}
		TEST_BEGIN("binary fixed point multiplication is correct") {
			ep2_rand(p);
			ep2_mul_pre_basic(t, p);
			bn_zero(k);
			ep2_mul_fix_basic(r, t, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_fix_basic(r, t, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep2_mul(r, p, k);
			ep2_mul_fix_basic(q, t, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_fix_basic(r, t, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE_BASIC; i++) {
			ep2_free(t[i]);
		}
#endif

#if EP_FIX == COMBS || !defined(STRIP)
		for (int i = 0; i < RLC_EP_TABLE_COMBS; i++) {
			ep2_new(t[i]);
		}
		TEST_BEGIN("single-table comb fixed point multiplication is correct") {
			ep2_rand(p);
			ep2_mul_pre_combs(t, p);
			bn_zero(k);
			ep2_mul_fix_combs(r, t, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_fix_combs(r, t, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep2_mul(r, p, k);
			ep2_mul_fix_combs(q, t, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_fix_combs(r, t, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE_COMBS; i++) {
			ep2_free(t[i]);
		}
#endif

#if EP_FIX == COMBD || !defined(STRIP)
		for (int i = 0; i < RLC_EP_TABLE_COMBD; i++) {
			ep2_new(t[i]);
		}
		TEST_BEGIN("double-table comb fixed point multiplication is correct") {
			ep2_rand(p);
			ep2_mul_pre_combd(t, p);
			bn_zero(k);
			ep2_mul_fix_combd(r, t, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_fix_combd(r, t, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep2_mul(r, p, k);
			ep2_mul_fix_combd(q, t, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_fix_combd(r, t, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE_COMBD; i++) {
			ep2_free(t[i]);
		}
#endif

#if EP_FIX == LWNAF || !defined(STRIP)
		for (int i = 0; i < RLC_EP_TABLE_LWNAF; i++) {
			ep2_new(t[i]);
		}
		TEST_BEGIN("left-to-right w-naf fixed point multiplication is correct") {
			ep2_rand(p);
			ep2_mul_pre_lwnaf(t, p);
			bn_zero(k);
			ep2_mul_fix_lwnaf(r, t, k);
			TEST_ASSERT(ep2_is_infty(r), end);
			bn_set_dig(k, 1);
			ep2_mul_fix_lwnaf(r, t, k);
			TEST_ASSERT(ep2_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep2_mul(r, p, k);
			ep2_mul_fix_lwnaf(q, t, k);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_fix_lwnaf(r, t, k);
			ep2_neg(r, r);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE_LWNAF; i++) {
			ep2_free(t[i]);
		}
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	ep2_free(p);
	ep2_free(q);
	ep2_free(r);
	bn_free(n);
	bn_free(k);
	return code;
}

static int simultaneous(void) {
	int code = RLC_ERR;
	bn_t n, k, l;
	ep2_t p, q, r;

	bn_null(n);
	bn_null(k);
	bn_null(l);
	ep2_null(p);
	ep2_null(q);
	ep2_null(r);

	TRY {
		bn_new(n);
		bn_new(k);
		bn_new(l);
		ep2_new(p);
		ep2_new(q);
		ep2_new(r);

		ep2_curve_get_gen(p);
		ep2_curve_get_ord(n);

		TEST_BEGIN("simultaneous point multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep2_mul(q, p, l);
			ep2_mul_sim(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep2_mul(q, p, k);
			ep2_mul_sim(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep2_mul_sim(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_sim(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep2_mul_sim(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;

#if EP_SIM == BASIC || !defined(STRIP)
		TEST_BEGIN("basic simultaneous point multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep2_mul(q, p, l);
			ep2_mul_sim_basic(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep2_mul(q, p, k);
			ep2_mul_sim_basic(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep2_mul_sim_basic(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_sim_basic(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep2_mul_sim_basic(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_SIM == TRICK || !defined(STRIP)
		TEST_BEGIN("shamir's trick for simultaneous multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep2_mul(q, p, l);
			ep2_mul_sim_trick(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep2_mul(q, p, k);
			ep2_mul_sim_trick(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep2_mul_sim_trick(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_sim_trick(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep2_mul_sim_trick(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_SIM == INTER || !defined(STRIP)
		TEST_BEGIN("interleaving for simultaneous multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep2_mul(q, p, l);
			ep2_mul_sim_inter(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep2_mul(q, p, k);
			ep2_mul_sim_inter(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep2_mul_sim_inter(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_sim_inter(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep2_mul_sim_inter(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_SIM == JOINT || !defined(STRIP)
		TEST_BEGIN("jsf for simultaneous multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep2_mul(q, p, l);
			ep2_mul_sim_joint(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep2_mul(q, p, k);
			ep2_mul_sim_joint(r, p, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep2_mul_sim_joint(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_sim_joint(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep2_mul_sim_joint(r, p, k, q, l);
			ep2_mul(p, p, k);
			ep2_mul(q, q, l);
			ep2_add(q, q, p);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("simultaneous multiplication with generator is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep2_mul(q, p, l);
			ep2_mul_sim_gen(r, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep2_mul_gen(q, k);
			ep2_mul_sim_gen(r, k, p, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep2_mul_sim_gen(r, k, q, l);
			ep2_curve_get_gen(p);
			ep2_mul_sim(q, p, k, q, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep2_mul_sim_gen(r, k, q, l);
			ep2_curve_get_gen(p);
			ep2_mul_sim(q, p, k, q, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep2_mul_sim_gen(r, k, q, l);
			ep2_curve_get_gen(p);
			ep2_mul_sim(q, p, k, q, l);
			TEST_ASSERT(ep2_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	bn_free(l);
	ep2_free(p);
	ep2_free(q);
	ep2_free(r);
	return code;
}

static int compression(void) {
	int code = RLC_ERR;
	ep2_t a, b, c;

	ep2_null(a);
	ep2_null(b);
	ep2_null(c);

	TRY {
		ep2_new(a);
		ep2_new(b);
		ep2_new(c);

		TEST_BEGIN("point compression is correct") {
			ep2_rand(a);
			ep2_pck(b, a);
			TEST_ASSERT(ep2_upk(c, b) == 1, end);
			TEST_ASSERT(ep2_cmp(a, c) == RLC_EQ, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		ERROR(end);
	}
	code = RLC_OK;
  end:
	ep2_free(a);
	ep2_free(b);
	ep2_free(c);
	return code;
}

static int hashing(void) {
	int code = RLC_ERR;
	bn_t n;
	ep2_t p;
	uint8_t msg[5];

	bn_null(n);
	ep2_null(p);

	TRY {
		bn_new(n);
		ep2_new(p);

		ep2_curve_get_ord(n);

		TEST_BEGIN("point hashing is correct") {
			rand_bytes(msg, sizeof(msg));
			ep2_map(p, msg, sizeof(msg));
			ep2_mul(p, p, n);
			TEST_ASSERT(ep2_is_infty(p) == 1, end);
		}
		TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	ep2_free(p);
	return code;
}

static int frobenius(void) {
	int code = RLC_ERR;
	ep2_t a, b, c;
	bn_t d, n;

	ep2_null(a);
	ep2_null(b);
	ep2_null(c);
	bn_null(d);
	bn_null(n);

	TRY {
		ep2_new(a);
		ep2_new(b);
		ep2_new(c);
		bn_new(d);
		bn_new(n);

		ep2_curve_get_ord(n);

		TEST_BEGIN("frobenius and scalar multiplication are consistent") {
			ep2_rand(a);
			ep2_frb(b, a, 1);
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			bn_mod(d, d, n);
			ep2_mul(c, a, d);
			TEST_ASSERT(ep2_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("squared frobenius is consistent") {
			ep2_rand(a);
			ep2_frb(b, a, 1);
			ep2_frb(b, b, 1);
			ep2_frb(c, a, 2);
			TEST_ASSERT(ep2_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("cubed frobenius is consistent") {
			ep2_rand(a);
			ep2_frb(b, a, 1);
			ep2_frb(b, b, 1);
			ep2_frb(b, b, 1);
			ep2_frb(c, a, 3);
			TEST_ASSERT(ep2_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	ep2_free(a);
	ep2_free(b);
	ep2_free(c);
	bn_free(d);
	bn_free(n);
	return code;
}
int util3(void) {
	int l, code = RLC_ERR;
	ep3_t a, b, c;
	bn_t n;
	uint8_t bin[4 * RLC_FP_BYTES + 1];

	ep3_null(a);
	ep3_null(b);
	ep3_null(c);
	bn_null(n);

	TRY{
		ep3_new(a);
		ep3_new(b);
		ep3_new(c);
		bn_new(n);
		TEST_BEGIN("comparison is consistent") {
			ep3_rand(a);
			ep3_rand(b);
			TEST_ASSERT(ep3_cmp(a, b) != RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("copy and comparison are consistent") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_rand(c);
			/* Compare points in affine coordinates. */
			if (ep3_cmp(a, c) != RLC_EQ) {
				ep3_copy(c, a);
				TEST_ASSERT(ep3_cmp(c, a) == RLC_EQ, end);
			}
			if (ep3_cmp(b, c) != RLC_EQ) {
				ep3_copy(c, b);
				TEST_ASSERT(ep3_cmp(b, c) == RLC_EQ, end);
			}
			/* Compare with one point in projective. */
			ep3_dbl(c, a);
			ep3_norm(c, c);
			ep3_dbl(a, a);
			TEST_ASSERT(ep3_cmp(c, a) == RLC_EQ, end);
			TEST_ASSERT(ep3_cmp(a, c) == RLC_EQ, end);
			/* Compare with two points in projective. */
			ep3_dbl(c, c);
			ep3_dbl(a, a);
			TEST_ASSERT(ep3_cmp(c, a) == RLC_EQ, end);
			TEST_ASSERT(ep3_cmp(a, c) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("negation and comparison are consistent") {
			ep3_rand(a);
			ep3_neg(b, a);
			TEST_ASSERT(ep3_cmp(a, b) != RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN
				("assignment to random/infinity and comparison are consistent")
		{
			ep3_rand(a);
			ep3_set_infty(c);
			TEST_ASSERT(ep3_cmp(a, c) != RLC_EQ, end);
			TEST_ASSERT(ep3_cmp(c, a) != RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("assignment to infinity and infinity test are consistent") {
			ep3_set_infty(a);
			TEST_ASSERT(ep3_is_infty(a), end);
		}
		TEST_END;
		TEST_BEGIN("validity test is correct") {
			ep3_set_infty(a);
			TEST_ASSERT(ep3_is_valid(a), end);
			ep3_rand(a);
			TEST_ASSERT(ep3_is_valid(a), end);
			fp3_rand(a->x);
			TEST_ASSERT(!ep3_is_valid(a), end);
		}
		TEST_END;

		/*TEST_BEGIN("reading and writing a point are consistent") {
			for (int j = 0; j < 2; j++) {
				ep3_set_infty(a);
				l = ep3_size_bin(a, j);
				ep3_write_bin(bin, l, a, j);
				ep3_read_bin(b, bin, l);
				TEST_ASSERT(ep3_cmp(a, b) == RLC_EQ, end);
				ep3_rand(a);
				l = ep3_size_bin(a, j);
				ep3_write_bin(bin, l, a, j);
				ep3_read_bin(b, bin, l);
				TEST_ASSERT(ep3_cmp(a, b) == RLC_EQ, end);
				ep3_rand(a);
				ep3_dbl(a, a);
				l = ep3_size_bin(a, j);
				ep3_norm(a, a);
				ep3_write_bin(bin, l, a, j);
				ep3_read_bin(b, bin, l);
				TEST_ASSERT(ep3_cmp(a, b) == RLC_EQ, end);
			}
		}
		TEST_END;*/
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	ep3_free(a);
	ep3_free(b);
	ep3_free(c);
	bn_free(n);
	return code;
}

int addition3(void) {
	int code = RLC_ERR;
	ep3_t a, b, c, d, e;

	ep3_null(a);
	ep3_null(b);
	ep3_null(c);
	ep3_null(d);
	ep3_null(e);

	TRY{
		ep3_new(a);
		ep3_new(b);
		ep3_new(c);
		ep3_new(d);
		ep3_new(e);

		TEST_BEGIN("point addition is commutative") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_add(d, a, b);
			ep3_add(e, b, a);
			TEST_ASSERT(ep3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point addition is associative") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_rand(c);
			ep3_add(d, a, b);
			ep3_add(d, d, c);
			ep3_add(e, b, c);
			ep3_add(e, e, a);
			TEST_ASSERT(ep3_cmp(d, e) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point addition has identity") {
			ep3_rand(a);
			ep3_set_infty(d);
			ep3_add(e, a, d);
			TEST_ASSERT(ep3_cmp(e, a) == RLC_EQ, end);
			ep3_add(e, d, a);
			TEST_ASSERT(ep3_cmp(e, a) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point addition has inverse") {
			ep3_rand(a);
			ep3_neg(d, a);
			ep3_add(e, a, d);
			TEST_ASSERT(ep3_is_infty(e), end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("point addition in affine coordinates is correct") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_add(d, a, b);
			ep3_norm(d, d);
			ep3_add_basic(e, a, b);
			TEST_ASSERT(ep3_cmp(e, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
#if !defined(EP_MIXED) || !defined(STRIP)
		TEST_BEGIN("point addition in projective coordinates is correct") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_add_projc(a, a, b);
			ep3_rand(b);
			ep3_rand(c);
			ep3_add_projc(b, b, c);
			ep3_norm(b, b);
			/* a and b in projective coordinates. */
			ep3_add_projc(d, a, b);
			ep3_norm(d, d);
			ep3_norm(a, a);
			ep3_norm(b, b);
			ep3_add(e, a, b);
			ep3_norm(e, e);
			TEST_ASSERT(ep3_cmp(e, d) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("point addition in mixed coordinates (z2 = 1) is correct") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_add_projc(a, a, b);
			ep3_rand(b);
			/* a and b in projective coordinates. */
			ep3_add_projc(d, a, b);
			ep3_norm(d, d);
			/* a in affine coordinates. */
			ep3_norm(a, a);
			ep3_add(e, a, b);
			ep3_norm(e, e);
			TEST_ASSERT(ep3_cmp(e, d) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point addition in mixed coordinates (z1,z2 = 1) is correct") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_norm(a, a);
			ep3_norm(b, b);
			/* a and b in affine coordinates. */
			ep3_add(d, a, b);
			ep3_norm(d, d);
			ep3_add_projc(e, a, b);
			ep3_norm(e, e);
			TEST_ASSERT(ep3_cmp(e, d) == RLC_EQ, end);
		} TEST_END;
#endif

	}
		CATCH_ANY{
			ERROR(end);
	}
	code = RLC_OK;
end:
	ep3_free(a);
	ep3_free(b);
	ep3_free(c);
	ep3_free(d);
	ep3_free(e);
	return code;
}

int subtraction3(void) {
	int code = RLC_ERR;
	ep3_t a, b, c, d;

	ep3_null(a);
	ep3_null(b);
	ep3_null(c);
	ep3_null(d);

	TRY{
		ep3_new(a);
		ep3_new(b);
		ep3_new(c);
		ep3_new(d);

		TEST_BEGIN("point subtraction is anti-commutative") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_sub(c, a, b);
			ep3_sub(d, b, a);
			ep3_neg(d, d);
			TEST_ASSERT(ep3_cmp(c, d) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("point subtraction has identity") {
			ep3_rand(a);
			ep3_set_infty(c);
			ep3_sub(d, a, c);
			TEST_ASSERT(ep3_cmp(d, a) == RLC_EQ, end);
		}
		TEST_END;

		TEST_BEGIN("point subtraction has inverse") {
			ep3_rand(a);
			ep3_sub(c, a, a);
			TEST_ASSERT(ep3_is_infty(c), end);
		}
		TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("point subtraction in affine coordinates is correct") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_sub(c, a, b);
			ep3_norm(c, c);
			ep3_sub_basic(d, a, b);
			TEST_ASSERT(ep3_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
#if !defined(EP_MIXED) || !defined(STRIP)
		TEST_BEGIN("point subtraction in projective coordinates is correct") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_add_projc(a, a, b);
			ep3_rand(b);
			ep3_rand(c);
			ep3_add_projc(b, b, c);
			/* a and b in projective coordinates. */
			ep3_sub_projc(c, a, b);
			ep3_norm(c, c);
			ep3_norm(a, a);
			ep3_norm(b, b);
			ep3_sub(d, a, b);
			ep3_norm(d, d);
			TEST_ASSERT(ep3_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("point subtraction in mixed coordinates (z2 = 1) is correct") {
			ep3_rand(a);
			ep3_rand(b);
			ep3_add_projc(a, a, b);
			ep3_rand(b);
			/* a and b in projective coordinates. */
			ep3_sub_projc(c, a, b);
			ep3_norm(c, c);
			/* a in affine coordinates. */
			ep3_norm(a, a);
			ep3_sub(d, a, b);
			ep3_norm(d, d);
			TEST_ASSERT(ep3_cmp(c, d) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN
				("point subtraction in mixed coordinates (z1,z2 = 1) is correct")
		{
			ep3_rand(a);
			ep3_rand(b);
			ep3_norm(a, a);
			ep3_norm(b, b);
			/* a and b in affine coordinates. */
			ep3_sub(c, a, b);
			ep3_norm(c, c);
			ep3_sub_projc(d, a, b);
			ep3_norm(d, d);
			TEST_ASSERT(ep3_cmp(c, d) == RLC_EQ, end);
		} TEST_END;
#endif
	}
		CATCH_ANY{
			ERROR(end);
	}
	code = RLC_OK;
end:
	ep3_free(a);
	ep3_free(b);
	ep3_free(c);
	ep3_free(d);
	return code;
}

int doubling3(void) {
	int code = RLC_ERR;
	ep3_t a, b, c;

	ep3_null(a);
	ep3_null(b);
	ep3_null(c);

	TRY{
		ep3_new(a);
		ep3_new(b);
		ep3_new(c);

		TEST_BEGIN("point doubling is correct") {
			ep3_rand(a);
			ep3_add(b, a, a);
			ep3_dbl(c, a);
			TEST_ASSERT(ep3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("point doubling in affine coordinates is correct") {
			ep3_rand(a);
			ep3_dbl(b, a);
			ep3_norm(b, b);
			ep3_dbl_basic(c, a);
			TEST_ASSERT(ep3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
		TEST_BEGIN("point doubling in projective coordinates is correct") {
			ep3_rand(a);
			ep3_dbl_projc(a, a);
			/* a in projective coordinates. */
			ep3_dbl_projc(b, a);
			ep3_norm(b, b);
			ep3_norm(a, a);
			ep3_dbl(c, a);
			ep3_norm(c, c);
			TEST_ASSERT(ep3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("point doubling in mixed coordinates (z1 = 1) is correct") {
			ep3_rand(a);
			ep3_dbl_projc(b, a);
			ep3_norm(b, b);
			ep3_dbl(c, a);
			ep3_norm(c, c);
			TEST_ASSERT(ep3_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
#endif
	}
		CATCH_ANY{
			ERROR(end);
	}
	code = RLC_OK;
end:
	ep3_free(a);
	ep3_free(b);
	ep3_free(c);
	return code;
}

static int multiplication3(void) {
	int code = RLC_ERR;
	bn_t n, k;
	ep3_t p, q, r;

	bn_null(n);
	bn_null(k);
	ep3_null(p);
	ep3_null(q);
	ep3_null(r);

	TRY{
		bn_new(n);
		bn_new(k);
		ep3_new(p);
		ep3_new(q);
		ep3_new(r);

		ep3_curve_get_gen(p);
		ep3_curve_get_ord(n);
		TEST_BEGIN("generator has the right order") {
			TEST_ASSERT(ep3_is_valid(p), end);
			ep3_mul(r, p, n);
			TEST_ASSERT(ep3_is_infty(r) == 1, end);
		} TEST_END;

		TEST_BEGIN("generator multiplication is correct") {
			bn_zero(k);
			ep3_mul_gen(r, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_gen(r, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep3_mul_basic(q, p, k);
			ep3_mul_gen(r, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_gen(r, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;

#if EP_MUL == BASIC || !defined(STRIP)
		TEST_BEGIN("binary point multiplication is correct") {
			bn_zero(k);
			ep3_mul_basic(r, p, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_basic(r, p, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			ep3_rand(p);
			ep3_mul_basic(r, p, n);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_rand_mod(k, n);
			ep3_mul_basic(q, p, k);
			ep3_mul_basic(r, p, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_basic(r, p, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_MUL == MONTY || !defined(STRIP)
		TEST_BEGIN("sliding window point multiplication is correct") {
			bn_zero(k);
			ep3_mul_slide(r, p, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_slide(r, p, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			ep3_rand(p);
			ep3_mul_basic(r, p, n);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_rand_mod(k, n);
			ep3_mul_basic(q, p, k);
			ep3_mul_slide(r, p, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_slide(r, p, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		}
		TEST_END;
#endif

#if EP_MUL == MONTY || !defined(STRIP)
		TEST_BEGIN("montgomery laddering point multiplication is correct") {
			bn_zero(k);
			ep3_mul_monty(r, p, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_monty(r, p, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			ep3_rand(p);
			ep3_mul_basic(r, p, n);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_rand_mod(k, n);
			ep3_mul_basic(q, p, k);
			ep3_mul_monty(r, p, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_monty(r, p, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		}
		TEST_END;
#endif

#if EP_MUL == LWNAF || !defined(STRIP)
		/*TEST_BEGIN("left-to-right w-naf point multiplication is correct") {
			bn_zero(k);
			ep3_mul_lwnaf(r, p, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_lwnaf(r, p, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			ep3_rand(p);
			ep3_mul_basic(r, p, n);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_rand_mod(k, n);
			ep3_mul_basic(q, p, k);
			ep3_mul_lwnaf(r, p, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_lwnaf(r, p, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		}
		TEST_END;*/
#endif

		TEST_BEGIN("multiplication by digit is correct") {
			ep3_mul_dig(r, p, 0);
			TEST_ASSERT(ep3_is_infty(r), end);
			ep3_mul_dig(r, p, 1);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			bn_rand(k, RLC_POS, RLC_DIG);
			ep3_mul_basic(q, p, k);
			ep3_mul_dig(r, p, k->dp[0]);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		}
		TEST_END;
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	bn_free(n);
	bn_free(k);
	ep3_free(p);
	ep3_free(q);
	ep3_free(r);
	return code;
}

static int fixed3(void) {
	int code = RLC_ERR;
	bn_t n, k;
	ep3_t p, q, r, t[RLC_EPX_TABLE_MAX];

	bn_null(n);
	bn_null(k);
	ep3_null(p);
	ep3_null(q);
	ep3_null(r);

	for (int i = 0; i < RLC_EPX_TABLE_MAX; i++) {
		ep3_null(t[i]);
	}

	TRY{
		bn_new(n);
		bn_new(k);
		ep3_new(p);
		ep3_new(q);
		ep3_new(r);

		ep3_curve_get_gen(p);
		ep3_curve_get_ord(n);

		for (int i = 0; i < RLC_EP_TABLE; i++) {
			ep3_new(t[i]);
		}
		TEST_BEGIN("fixed point multiplication is correct") {
			ep3_rand(p);
			ep3_mul_pre(t, p);
			bn_zero(k);
			ep3_mul_fix(r, t, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_fix(r, t, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep3_mul(q, p, k);
			ep3_mul_fix(q, t, k);
			ep3_mul(r, p, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_fix(r, t, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE; i++) {
			ep3_free(t[i]);
		}

#if EP_FIX == BASIC || !defined(STRIP)
		for (int i = 0; i < RLC_EP_TABLE_BASIC; i++) {
			ep3_new(t[i]);
		}
		TEST_BEGIN("binary fixed point multiplication is correct") {
			ep3_rand(p);
			ep3_mul_pre_basic(t, p);
			bn_zero(k);
			ep3_mul_fix_basic(r, t, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_fix_basic(r, t, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep3_mul(r, p, k);
			ep3_mul_fix_basic(q, t, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_fix_basic(r, t, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE_BASIC; i++) {
			ep3_free(t[i]);
		}
#endif

#if EP_FIX == COMBS || !defined(STRIP)
		for (int i = 0; i < RLC_EP_TABLE_COMBS; i++) {
			ep3_new(t[i]);
		}
		TEST_BEGIN("single-table comb fixed point multiplication is correct") {
			ep3_rand(p);
			ep3_mul_pre_combs(t, p);
			bn_zero(k);
			ep3_mul_fix_combs(r, t, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_fix_combs(r, t, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep3_mul(r, p, k);
			ep3_mul_fix_combs(q, t, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_fix_combs(r, t, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE_COMBS; i++) {
			ep3_free(t[i]);
		}
#endif

#if EP_FIX == COMBD || !defined(STRIP)
		for (int i = 0; i < RLC_EP_TABLE_COMBD; i++) {
			ep3_new(t[i]);
		}
		TEST_BEGIN("double-table comb fixed point multiplication is correct") {
			ep3_rand(p);
			ep3_mul_pre_combd(t, p);
			bn_zero(k);
			ep3_mul_fix_combd(r, t, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_fix_combd(r, t, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep3_mul(r, p, k);
			ep3_mul_fix_combd(q, t, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_fix_combd(r, t, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE_COMBD; i++) {
			ep3_free(t[i]);
		}
#endif

#if EP_FIX == LWNAF || !defined(STRIP)
		for (int i = 0; i < RLC_EP_TABLE_LWNAF; i++) {
			ep3_new(t[i]);
		}
		TEST_BEGIN("left-to-right w-naf fixed point multiplication is correct") {
			ep3_rand(p);
			ep3_mul_pre_lwnaf(t, p);
			bn_zero(k);
			ep3_mul_fix_lwnaf(r, t, k);
			TEST_ASSERT(ep3_is_infty(r), end);
			bn_set_dig(k, 1);
			ep3_mul_fix_lwnaf(r, t, k);
			TEST_ASSERT(ep3_cmp(p, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			ep3_mul(r, p, k);
			ep3_mul_fix_lwnaf(q, t, k);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_fix_lwnaf(r, t, k);
			ep3_neg(r, r);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
		for (int i = 0; i < RLC_EP_TABLE_LWNAF; i++) {
			ep3_free(t[i]);
		}
#endif
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	ep3_free(p);
	ep3_free(q);
	ep3_free(r);
	bn_free(n);
	bn_free(k);
	return code;
}

static int simultaneous3(void) {
	int code = RLC_ERR;
	bn_t n, k, l;
	ep3_t p, q, r;

	bn_null(n);
	bn_null(k);
	bn_null(l);
	ep3_null(p);
	ep3_null(q);
	ep3_null(r);

	TRY{
		bn_new(n);
		bn_new(k);
		bn_new(l);
		ep3_new(p);
		ep3_new(q);
		ep3_new(r);

		ep3_curve_get_gen(p);
		ep3_curve_get_ord(n);

		TEST_BEGIN("simultaneous point multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep3_mul(q, p, l);
			ep3_mul_sim(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep3_mul(q, p, k);
			ep3_mul_sim(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep3_mul_sim(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_sim(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep3_mul_sim(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;

#if EP_SIM == BASIC || !defined(STRIP)
		TEST_BEGIN("basic simultaneous point multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep3_mul(q, p, l);
			ep3_mul_sim_basic(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep3_mul(q, p, k);
			ep3_mul_sim_basic(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep3_mul_sim_basic(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_sim_basic(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep3_mul_sim_basic(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_SIM == TRICK || !defined(STRIP)
		TEST_BEGIN("shamir's trick for simultaneous multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep3_mul(q, p, l);
			ep3_mul_sim_trick(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep3_mul(q, p, k);
			ep3_mul_sim_trick(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep3_mul_sim_trick(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_sim_trick(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep3_mul_sim_trick(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_SIM == INTER || !defined(STRIP)
		TEST_BEGIN("interleaving for simultaneous multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep3_mul(q, p, l);
			ep3_mul_sim_inter(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep3_mul(q, p, k);
			ep3_mul_sim_inter(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep3_mul_sim_inter(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_sim_inter(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep3_mul_sim_inter(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_SIM == JOINT || !defined(STRIP)
		TEST_BEGIN("jsf for simultaneous multiplication is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep3_mul(q, p, l);
			ep3_mul_sim_joint(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep3_mul(q, p, k);
			ep3_mul_sim_joint(r, p, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep3_mul_sim_joint(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_sim_joint(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep3_mul_sim_joint(r, p, k, q, l);
			ep3_mul(p, p, k);
			ep3_mul(q, q, l);
			ep3_add(q, q, p);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
#endif

		TEST_BEGIN("simultaneous multiplication with generator is correct") {
			bn_zero(k);
			bn_rand_mod(l, n);
			ep3_mul(q, p, l);
			ep3_mul_sim_gen(r, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_zero(l);
			ep3_mul_gen(q, k);
			ep3_mul_sim_gen(r, k, p, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_rand_mod(k, n);
			bn_rand_mod(l, n);
			ep3_mul_sim_gen(r, k, q, l);
			ep3_curve_get_gen(p);
			ep3_mul_sim(q, p, k, q, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(k, k);
			ep3_mul_sim_gen(r, k, q, l);
			ep3_curve_get_gen(p);
			ep3_mul_sim(q, p, k, q, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
			bn_neg(l, l);
			ep3_mul_sim_gen(r, k, q, l);
			ep3_curve_get_gen(p);
			ep3_mul_sim(q, p, k, q, l);
			TEST_ASSERT(ep3_cmp(q, r) == RLC_EQ, end);
		} TEST_END;
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	bn_free(n);
	bn_free(k);
	bn_free(l);
	ep3_free(p);
	ep3_free(q);
	ep3_free(r);
	return code;
}

static int compression3(void) {
	int code = RLC_ERR;
	ep3_t a, b, c;

	ep3_null(a);
	ep3_null(b);
	ep3_null(c);

	TRY{
		ep3_new(a);
		ep3_new(b);
		ep3_new(c);

		TEST_BEGIN("point compression is correct") {
			ep3_rand(a);
			ep3_pck(b, a);
			TEST_ASSERT(ep3_upk(c, b) == 1, end);
			TEST_ASSERT(ep3_cmp(a, c) == RLC_EQ, end);
		}
		TEST_END;
	}
		CATCH_ANY{
			ERROR(end);
	}
	code = RLC_OK;
end:
	ep3_free(a);
	ep3_free(b);
	ep3_free(c);
	return code;
}

static int hashing3(void) {
	int code = RLC_ERR;
	bn_t n;
	ep3_t p;
	uint8_t msg[5];

	bn_null(n);
	ep3_null(p);

	TRY{
		bn_new(n);
		ep3_new(p);

		ep3_curve_get_ord(n);

		TEST_BEGIN("point hashing is correct") {
			rand_bytes(msg, sizeof(msg));
			ep3_map(p, msg, sizeof(msg));
			ep3_mul(p, p, n);
			TEST_ASSERT(ep3_is_infty(p) == 1, end);
		}
		TEST_END;
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	bn_free(n);
	ep3_free(p);
	return code;
}

static int frobenius3(void) {
	int code = RLC_ERR;
	ep3_t a, b, c;
	bn_t d, n;

	ep3_null(a);
	ep3_null(b);
	ep3_null(c);
	bn_null(d);
	bn_null(n);

	TRY{
		ep3_new(a);
		ep3_new(b);
		ep3_new(c);
		bn_new(d);
		bn_new(n);

		ep3_curve_get_ord(n);

		TEST_BEGIN("frobenius and scalar multiplication are consistent") {
			ep3_rand(a);
			ep3_frb(b, a, 1);
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			bn_mod(d, d, n);
			ep3_mul_basic(c, a, d);
			TEST_ASSERT(ep3_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("squared frobenius is consistent") {
			ep3_rand(a);
			ep3_frb(b, a, 1);
			ep3_frb(b, b, 1);
			ep3_frb(c, a, 2);
			TEST_ASSERT(ep3_cmp(c, b) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("cubed frobenius is consistent") {
			ep3_rand(a);
			ep3_frb(b, a, 1);
			ep3_frb(b, b, 1);
			ep3_frb(b, b, 1);
			ep3_frb(c, a, 3);
			TEST_ASSERT(ep3_cmp(c, b) == RLC_EQ, end);
		} TEST_END;
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	ep3_free(a);
	ep3_free(b);
	ep3_free(c);
	bn_free(d);
	bn_free(n);
	return code;
}
int main(void) {
	int flag = 2;
	core_init();
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}
	util_banner("Tests for the EPX module", 0);
#if FP_PRIME==508 || FP_PRIME==676
	if (ep_param_set_any_endom() == RLC_ERR) {
		THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

#else
	if (ep_param_set_any_pairf() == RLC_ERR) {
		THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}
#endif
	ep_param_print();

#if FP_PRIME==508 || FP_PRIME==676
ep3_curve_set_twist(1);
	if (ep3_curve_is_twist() == 0) {
		THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}
	flag = 3;
#else
	if (ep2_curve_is_twist() == 0) {
		THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}
#endif

	util_banner("Sextic twist:", 0);
	util_banner("Utilities:", 1);

	if (flag == 2) {

		if (memory() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (util() != RLC_OK) {
			core_clean();
			return 1;
		}

		util_banner("Arithmetic:", 1);

		if (addition() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (subtraction() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (multiplication() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (fixed() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (simultaneous() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (compression() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (hashing() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (frobenius() != RLC_OK) {
			core_clean();
			return 1;
		}
	}
	else {
		if (memory() != RLC_OK) {
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

/*		if (fixed3() != RLC_OK) {
			core_clean();
			return 1;
		}*/

		if (simultaneous3() != RLC_OK) {
			core_clean();
			return 1;
		}

/*		if (compression3() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (hashing3() != RLC_OK) {
			core_clean();
			return 1;
		}*/

		if (frobenius3() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	util_banner("All tests have passed.\n", 0);

	core_clean();
	return 0;
}
