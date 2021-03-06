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
 * Implementation of point multiplication on prime elliptic curves over
 * cubic extensions.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if EP_MUL == LWNAF || !defined(STRIP)

#if defined(EP_ENDOM)
//question???
static void ep3_mul_glv_imp(ep3_t r, ep3_t p, const bn_t k) {
	int i, j, l;
	bn_t n, _k[4], u[4], v[4];
	ep3_t q[4];

	bn_null(n);

	TRY {
		bn_new(n);
		for (i = 0; i < 4; i++) {
			bn_null(u[i]);
			bn_null(v[i]);
			bn_null(_k[i]);
			ep3_null(q[i]);
			bn_new(u[i]);
			bn_new(v[i]);
			bn_new(_k[i]);
			ep3_new(q[i]);
		}

		ep3_curve_get_ord(n);

		switch (ep_curve_is_pairf()) {
			case EP_BN:
				ep3_curve_get_vs(v);

				for (i = 0; i < 4; i++) {
					bn_mul(v[i], v[i], k);
					bn_div(v[i], v[i], n);
					if (bn_sign(v[i]) == RLC_NEG) {
						bn_add_dig(v[i], v[i], 1);
					}
					bn_zero(_k[i]);
				}

				/* u0 = x + 1, u1 = 2x + 1, u2 = 2x, u3 = x - 1. */
				fp_prime_get_par(u[0]);
				bn_dbl(u[2], u[0]);
				bn_add_dig(u[1], u[2], 1);
				bn_sub_dig(u[3], u[0], 1);
				bn_add_dig(u[0], u[0], 1);
				bn_copy(_k[0], k);
				for (i = 0; i < 4; i++) {
					bn_mul(u[i], u[i], v[i]);
					bn_mod(u[i], u[i], n);
					bn_add(_k[0], _k[0], n);
					bn_sub(_k[0], _k[0], u[i]);
					bn_mod(_k[0], _k[0], n);
				}

				/* u0 = x, u1 = -x, u2 = 2x + 1, u3 = 4x + 2. */
				fp_prime_get_par(u[0]);
				bn_neg(u[1], u[0]);
				bn_dbl(u[2], u[0]);
				bn_add_dig(u[2], u[2], 1);
				bn_dbl(u[3], u[2]);
				for (i = 0; i < 4; i++) {
					bn_mul(u[i], u[i], v[i]);
					bn_mod(u[i], u[i], n);
					bn_add(_k[1], _k[1], n);
					bn_sub(_k[1], _k[1], u[i]);
					bn_mod(_k[1], _k[1], n);
				}

				/* u0 = x, u1 = -(x + 1), u2 = 2x + 1, u3 = -(2x - 1). */
				fp_prime_get_par(u[0]);
				bn_add_dig(u[1], u[0], 1);
				bn_neg(u[1], u[1]);
				bn_dbl(u[2], u[0]);
				bn_add_dig(u[2], u[2], 1);
				bn_sub_dig(u[3], u[2], 2);
				bn_neg(u[3], u[3]);
				for (i = 0; i < 4; i++) {
					bn_mul(u[i], u[i], v[i]);
					bn_mod(u[i], u[i], n);
					bn_add(_k[2], _k[2], n);
					bn_sub(_k[2], _k[2], u[i]);
					bn_mod(_k[2], _k[2], n);
				}

				/* u0 = -2x, u1 = -x, u2 = 2x + 1, u3 = x - 1. */
				fp_prime_get_par(u[1]);
				bn_dbl(u[0], u[1]);
				bn_neg(u[0], u[0]);
				bn_dbl(u[2], u[1]);
				bn_add_dig(u[2], u[2], 1);
				bn_sub_dig(u[3], u[1], 1);
				bn_neg(u[1], u[1]);
				for (i = 0; i < 4; i++) {
					bn_mul(u[i], u[i], v[i]);
					bn_mod(u[i], u[i], n);
					bn_add(_k[3], _k[3], n);
					bn_sub(_k[3], _k[3], u[i]);
					bn_mod(_k[3], _k[3], n);
				}

				for (i = 0; i < 4; i++) {
					l = bn_bits(_k[i]);
					bn_sub(_k[i], n, _k[i]);
					if (bn_bits(_k[i]) > l) {
						bn_sub(_k[i], _k[i], n);
						_k[i]->sign = RLC_POS;
					} else {
						_k[i]->sign = RLC_NEG;
					}
				}
				break;
			case EP_OT:
			case EP_B12:
				bn_abs(v[0], k);
				fp_prime_get_par(u[0]);
				bn_copy(u[1], u[0]);
				if (bn_sign(u[0]) == RLC_NEG) {
					bn_neg(u[0], u[0]);
				}

				for (i = 0; i < 4; i++) {
					bn_mod(_k[i], v[0], u[0]);
					bn_div(v[0], v[0], u[0]);
					if ((bn_sign(u[1]) == RLC_NEG) && (i % 2 != 0)) {
						bn_neg(_k[i], _k[i]);
					}
					if (bn_sign(k) == RLC_NEG) {
						bn_neg(_k[i], _k[i]);
					}
				}

				break;
		}

		ep3_norm(q[0], p);
		ep3_frb(q[1], q[0], 1);
		ep3_frb(q[2], q[1], 1);
		ep3_frb(q[3], q[2], 1);

		for (i = 0; i < 4; i++) {
			if (bn_sign(_k[i]) == RLC_NEG) {
				ep3_neg(q[i], q[i]);
			}
		}

		l = RLC_MAX(bn_bits(_k[0]), bn_bits(_k[1]));
		l = RLC_MAX(l, RLC_MAX(bn_bits(_k[2]), bn_bits(_k[3])));
		ep3_set_infty(r);
		for (i = l - 1; i >= 0; i--) {
			ep3_dbl(r, r);
			for (j = 0; j < 4; j++) {
				if (bn_get_bit(_k[j], i)) {
					ep3_add(r, r, q[j]);
				}
			}
		}

		/* Convert r to affine coordinates. */
		ep3_norm(r, r);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(n);
		for (i = 0; i < 4; i++) {
			bn_free(u[i]);
			bn_free(v[i]);
			bn_free(_k[i]);
			ep3_free(q[i]);
		}

	}
}

#endif /* EP_ENDOM */

#if defined(EP_PLAIN) || defined(EP_SUPER)

static void ep3_mul_naf_imp(ep3_t r, ep3_t p, const bn_t k) {
	int l, i, n;
	int8_t naf[RLC_FP_BITS + 1];
	ep3_t t[1 << (EP_WIDTH - 2)];

	TRY {
		/* Prepare the precomputation table. */
		for (i = 0; i < (1 << (EP_WIDTH - 2)); i++) {
			ep3_null(t[i]);
			ep3_new(t[i]);
		}
		/* Compute the precomputation table. */
		ep3_tab(t, p, EP_WIDTH);

		/* Compute the w-NAF representation of k. */
		l = sizeof(naf);
		bn_rec_naf(naf, &l, k, EP_WIDTH);

		ep3_set_infty(r);
		for (i = l - 1; i >= 0; i--) {
			ep3_dbl(r, r);

			n = naf[i];
			if (n > 0) {
				ep3_add(r, r, t[n / 2]);
			}
			if (n < 0) {
				ep3_sub(r, r, t[-n / 2]);
			}
		}
		/* Convert r to affine coordinates. */
		ep3_norm(r, r);
		if (bn_sign(k) == RLC_NEG) {
			ep3_neg(r, r);
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		/* Free the precomputation table. */
		for (i = 0; i < (1 << (EP_WIDTH - 2)); i++) {
			ep3_free(t[i]);
		}
	}
}

#endif /* EP_PLAIN || EP_SUPER */
#endif /* EP_MUL == LWNAF */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

#if EP_MUL == BASIC || EP_MUL == LWNAF || !defined(STRIP)

void ep3_mul_basic(ep3_t r, ep3_t p, const bn_t k) {
	int i, l;
	ep3_t t;

	ep3_null(t);

	if (bn_is_zero(k) || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	TRY {
		ep3_new(t);
		l = bn_bits(k);

		if (bn_get_bit(k, l - 1)) {
			ep3_copy(t, p);
		} else {
			ep3_set_infty(t);
		}

		for (i = l - 2; i >= 0; i--) {
			ep3_dbl(t, t);
			if (bn_get_bit(k, i)) {
				ep3_add(t, t, p);
			}
		}

		ep3_copy(r, t);
		ep3_norm(r, r);
		if (bn_sign(k) == RLC_NEG) {
			ep3_neg(r, r);
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep3_free(t);
	}
}

#endif

#if EP_MUL == SLIDE || !defined(STRIP)

void ep3_mul_slide(ep3_t r, ep3_t p, const bn_t k) {
	ep3_t t[1 << (EP_WIDTH - 1)], q;
	int i, j, l;
	uint8_t win[RLC_FP_BITS + 1];

	ep3_null(q);

	if (bn_is_zero(k) || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	TRY {
		for (i = 0; i < (1 << (EP_WIDTH - 1)); i ++) {
			ep3_null(t[i]);
			ep3_new(t[i]);
		}

		ep3_new(q);

		ep3_copy(t[0], p);
		ep3_dbl(q, p);

#if defined(EP_MIXED)
		ep3_norm(q, q);
#endif

		/* Create table. */
		for (i = 1; i < (1 << (EP_WIDTH - 1)); i++) {
			ep3_add(t[i], t[i - 1], q);
		}

#if defined(EP_MIXED)
		ep3_norm_sim(t + 1, t + 1, (1 << (EP_WIDTH - 1)) - 1);
#endif

		ep3_set_infty(q);
		l = RLC_FP_BITS + 1;
		bn_rec_slw(win, &l, k, EP_WIDTH);
		for (i = 0; i < l; i++) {
			if (win[i] == 0) {
				ep3_dbl(q, q);
			} else {
				for (j = 0; j < util_bits_dig(win[i]); j++) {
					ep3_dbl(q, q);
				}
				ep3_add(q, q, t[win[i] >> 1]);
			}
		}

		ep3_norm(r, q);
		if (bn_sign(k) == RLC_NEG) {
			ep3_neg(r, r);
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		for (i = 0; i < (1 << (EP_WIDTH - 1)); i++) {
			ep3_free(t[i]);
		}
		ep3_free(q);
	}
}

#endif

#if EP_MUL == MONTY || !defined(STRIP)
//may wrong
void ep3_mul_monty(ep3_t r, ep3_t p, const bn_t k) {
	ep3_t t[2];

	ep3_null(t[0]);
	ep3_null(t[1]);

	if (bn_is_zero(k) || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	TRY {
		ep3_new(t[0]);
		ep3_new(t[1]);

		ep3_set_infty(t[0]);
		ep3_copy(t[1], p);

		for (int i = bn_bits(k) - 1; i >= 0; i--) {
			int j = bn_get_bit(k, i);
			dv_swap_cond(t[0]->x[0], t[1]->x[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->x[1], t[1]->x[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->x[2], t[1]->x[2], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->y[0], t[1]->y[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->y[1], t[1]->y[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->y[2], t[1]->y[2], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->z[0], t[1]->z[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->z[1], t[1]->z[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->z[2], t[1]->z[2], RLC_FP_DIGS, j ^ 1);
			ep3_add(t[0], t[0], t[1]);
			ep3_dbl(t[1], t[1]);
			dv_swap_cond(t[0]->x[0], t[1]->x[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->x[1], t[1]->x[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->x[2], t[1]->x[2], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->y[0], t[1]->y[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->y[1], t[1]->y[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->y[2], t[1]->y[2], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->z[0], t[1]->z[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->z[1], t[1]->z[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->z[2], t[1]->z[2], RLC_FP_DIGS, j ^ 1);
		}

		ep3_norm(r, t[0]);
		if (bn_sign(k) == RLC_NEG) {
			ep3_neg(r, r);
		}
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep3_free(t[1]);
		ep3_free(t[0]);
	}
}

#endif

#if EP_MUL == LWNAF || !defined(STRIP)

void ep3_mul_lwnaf(ep3_t r, ep3_t p, const bn_t k) {
	if (bn_is_zero(k) || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

#if defined(EP_ENDOM)
	if (ep_curve_is_endom()) {
		ep3_mul_glv_imp(r, p, k);
		return;
	}
#endif

#if defined(EP_PLAIN) || defined(EP_SUPER)
	ep3_mul_naf_imp(r, p, k);
#endif
}

#endif

void ep3_mul_gen(ep3_t r, bn_t k) {
	if (bn_is_zero(k)) {
		ep3_set_infty(r);
		return;
	}

//#ifdef EP_PRECO
	//ep3_mul_fix(r, ep3_curve_get_tab(), k);
//#else
	ep3_t g;

	ep3_null(g);

	TRY {
		ep3_new(g);
		ep3_curve_get_gen(g);
		ep3_mul_basic(r, g, k);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep3_free(g);
	}
//#endif
}

void ep3_mul_dig(ep3_t r, ep3_t p, dig_t k) {
	int i, l;
	ep3_t t;

	ep3_null(t);

	if (k == 0 || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	TRY {
		ep3_new(t);

		l = util_bits_dig(k);

		ep3_copy(t, p);

		for (i = l - 2; i >= 0; i--) {
			ep3_dbl(t, t);
			if (k & ((dig_t)1 << i)) {
				ep3_add(t, t, p);
			}
		}

		ep3_norm(r, t);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep3_free(t);
	}
}
