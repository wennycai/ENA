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
 * Implementation of the pairings over prime curves.
 *
 * @ingroup pp
 */

#include "relic_core.h"
#include "relic_pp.h"
#include "relic_util.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Compute the Miller loop for pairings of type G_2 x G_1 over the bits of a
 * given parameter.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] p				- the first pairing argument in affine coordinates.
 * @param[in] q				- the second pairing argument in affine coordinates.
 * @param[in] m 			- the number of pairings to evaluate.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_k2(fp2_t r, ep_t *t, ep_t *p, ep_t *q, int m, bn_t a) {
	fp2_t l;
	ep_t *_q = RLC_ALLOCA(ep_t, m);
	int i, j;

	fp2_null(l);

	TRY {
		fp2_new(l);
		for (j = 0; j < m; j++) {
			ep_null(_q[j]);
			ep_new(_q[j]);
			ep_copy(t[j], p[j]);
			ep_neg(_q[j], q[j]);
		}

		fp2_zero(l);

		for (i = bn_bits(a) - 2; i >= 0; i--) {
			fp2_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_k2(l, t[j], t[j], _q[j]);
				fp2_mul(r, r, l);
				if (bn_get_bit(a, i)) {
					pp_add_k2(l, t[j], p[j], q[j]);
					fp2_mul(r, r, l);
				}
			}
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp2_free(l);
		for (j = 0; j < m; j++) {
			ep_free(_q[j]);
		}
	}
}

/**
 * Compute the Miller loop for pairings of type G_1 x G_2 over the bits of a
 * given parameter.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] p				- the first pairing argument in affine coordinates.
 * @param[in] q				- the second pairing argument in affine coordinates.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_lit_k2(fp2_t r, ep_t *t, ep_t *p, ep_t *q, int m, bn_t a) {
	fp2_t l, _l;
	ep_t *_q = RLC_ALLOCA(ep_t, m);
	int i, j;

	fp2_null(l);
	fp2_null(_l);

	TRY {
		fp2_new(l);
		fp2_new(_l);
		for (j = 0; j < m; j++) {
			ep_null(_q[j]);
			ep_new(_q[j]);
			ep_copy(t[j], p[j]);
			ep_neg(_q[j], q[j]);
		}

		for (i = bn_bits(a) - 2; i >= 0; i--) {
			fp2_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_k2(l, t[j], t[j], _q[j]);
				fp_copy(_l[0], l[1]);
				fp_copy(_l[1], l[0]);
				fp2_mul(r, r, _l);
				if (bn_get_bit(a, i)) {
					pp_add_k2(l, t[j], p[j], q[j]);
					fp_copy(_l[0], l[1]);
					fp_copy(_l[1], l[0]);
					fp2_mul(r, r, _l);
				}
			}
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp2_free(l);
		fp2_free(_l);
		for (j = 0; j < m; j++) {
			ep_null(_q[j]);
		}
	}
}

/**
 * Compute the Miller loop for pairings of type G_2 x G_1 over the bits of a
 * given parameter represented in sparse form.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] q				- the vector of first arguments in affine coordinates.
 * @param[in] p				- the vector of second arguments in affine coordinates.
 * @param[in] n 			- the number of pairings to evaluate.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_k12(fp12_t r, ep2_t *t, ep2_t *q, ep_t *p, int m, bn_t a) {
	fp12_t l;
	ep_t *_p = RLC_ALLOCA(ep_t, m);
	ep2_t *_q = RLC_ALLOCA(ep2_t, m);
	int i, j, len = bn_bits(a) + 1;
	int8_t *s = RLC_ALLOCA(int8_t, len);

	if (m == 0) {
		return;
	}

	fp12_null(l);

	TRY {
		fp12_new(l);
		fp12_zero(l);

		for (j = 0; j < m; j++) {
			ep_null(_p[j]);
			ep2_null(_q[j]);
			ep_new(_p[j]);
			ep2_new(_q[j]);
			ep2_copy(t[j], q[j]);
			ep2_neg(_q[j], q[j]);
#if EP_ADD == BASIC
			ep_neg(_p[j], p[j]);
#else
			fp_add(_p[j]->x, p[j]->x, p[j]->x);
			fp_add(_p[j]->x, _p[j]->x, p[j]->x);
			fp_neg(_p[j]->y, p[j]->y);
#endif
		}

		bn_rec_naf(s, &len, a, 2);
		pp_dbl_k12(r, t[0], t[0], _p[0]);
		for (j = 1; j < m; j++) {
			pp_dbl_k12(l, t[j], t[j], _p[j]);
			fp12_mul_dxs(r, r, l);
		}
		if (s[len - 2] > 0) {
			for (j = 0; j < m; j++) {
				pp_add_k12(l, t[j], q[j], p[j]);
				fp12_mul_dxs(r, r, l);
			}
		}
		if (s[len - 2] < 0) {
			for (j = 0; j < m; j++) {
				pp_add_k12(l, t[j], _q[j], p[j]);
				fp12_mul_dxs(r, r, l);
			}
		}

		for (i = len - 3; i >= 0; i--) {
			fp12_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_k12(l, t[j], t[j], _p[j]);
				fp12_mul_dxs(r, r, l);
				if (s[i] > 0) {
					pp_add_k12(l, t[j], q[j], p[j]);
					fp12_mul_dxs(r, r, l);
				}
				if (s[i] < 0) {
					pp_add_k12(l, t[j], _q[j], p[j]);
					fp12_mul_dxs(r, r, l);
				}
			}
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp12_free(l);
		for (j = 0; j < m; j++) {
			ep_free(_p[j]);
			ep2_free(_q[j]);
		}
	}
}

/**
 * Compute the Miller loop for pairings of type G_1 x G_2 over the bits of a
 * given parameter.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] p				- the first pairing argument in affine coordinates.
 * @param[in] q				- the second pairing argument in affine coordinates.
 * @param[in] n 			- the number of pairings to evaluate.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_lit_k12(fp12_t r, ep_t *t, ep_t *p, ep2_t *q, int m, bn_t a) {
	fp12_t l;
	ep2_t *_q = RLC_ALLOCA(ep2_t, m);
	int j;

	fp12_null(l);

	TRY {
		fp12_new(l);
		for (j = 0; j < m; j++) {
			ep2_null(_q[j]);
			ep2_new(_q[j]);
			ep_copy(t[j], p[j]);
			ep2_neg(_q[j], q[j]);
		}

		fp12_zero(l);
		for (int i = bn_bits(a) - 2; i >= 0; i--) {
			fp12_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_lit_k12(l, t[j], t[j], _q[j]);
				fp12_mul(r, r, l);
				if (bn_get_bit(a, i)) {
					pp_add_lit_k12(l, t[j], p[j], q[j]);
					fp12_mul(r, r, l);
				}
			}
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp12_free(l);
		for (j = 0; j < m; j++) {
			ep2_free(_q[j]);
		}
	}
}

/**
 * Compute the final lines for optimal ate pairings.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] q				- the first point of the pairing, in G_2.
 * @param[in] p				- the second point of the pairing, in G_1.
 * @param[in] a				- the loop parameter.
 */
static void pp_fin_k12_oatep(fp12_t r, ep2_t t, ep2_t q, ep_t p) {
	ep2_t q1, q2;
	fp12_t tmp;

	fp12_null(tmp);
	ep2_null(q1);
	ep2_null(q2);

	TRY {
		ep2_new(q1);
		ep2_new(q2);
		fp12_new(tmp);
		fp12_zero(tmp);

		fp2_set_dig(q1->z, 1);
		fp2_set_dig(q2->z, 1);

		ep2_frb(q1, q, 1);
		ep2_frb(q2, q, 2);
		ep2_neg(q2, q2);

		pp_add_k12(tmp, t, q1, p);
		fp12_mul_dxs(r, r, tmp);
		pp_add_k12(tmp, t, q2, p);
		fp12_mul_dxs(r, r, tmp);
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		fp12_free(tmp);
		ep2_free(q1);
		ep2_free(q2);
	}
}


/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void pp_map_init(void) {
	ep2_curve_init();
}

void pp_map_clean(void) {
	ep2_curve_clean();
}

#if PP_MAP == TATEP || PP_MAP == OATEP || !defined(STRIP)

void pp_map_tatep_k2(fp2_t r, ep_t p, ep_t q) {
	ep_t _p[1], _q[1], t[1];
	bn_t n;

	ep_null(_p[0]);
	ep_null(_q[0]);
	ep_null(t[0]);
	bn_null(n);

	TRY {
		ep_new(t[0]);
		bn_new(n);

		ep_norm(_p[0], p);
		ep_norm(_q[0], q);
		ep_curve_get_ord(n);
		/* Since p has order n, we do not have to perform last iteration. */
		bn_sub_dig(n, n, 1);
		fp2_set_dig(r, 1);

		if (!ep_is_infty(p) && !ep_is_infty(q)) {
			pp_mil_k2(r, t, _p, _q, 1, n);
			pp_exp_k2(r, r);
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep_free(_p[0]);
		ep_free(_q[0]);
		ep_free(t[0]);
		bn_free(n);
	}
}

void pp_map_sim_tatep_k2(fp2_t r, ep_t *p, ep_t *q, int m) {
	ep_t *_p =RLC_ALLOCA(ep_t, m),
        *_q = RLC_ALLOCA(ep_t, m),
        *t = RLC_ALLOCA(ep_t, m);
	bn_t n;
	int i, j;

	bn_null(n);

	TRY {
		bn_new(n);
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(_q[i]);
			ep_null(t[i]);
			ep_new(_p[i]);
			ep_new(_q[i]);
			ep_new(t[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep_norm(_q[j++], q[i]);
			}
		}

		ep_curve_get_ord(n);
		bn_sub_dig(n, n, 1);
		fp2_set_dig(r, 1);
		if (j > 0) {
			pp_mil_k2(r, t, _p, _q, j, n);
			pp_exp_k2(r, r);
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(n);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep_free(_q[i]);
			ep_free(t[i]);
		}
	}
}

#endif

#if PP_MAP == TATEP || !defined(STRIP)

void pp_map_tatep_k12(fp12_t r, ep_t p, ep2_t q) {
	ep_t _p[1], t[1];
	ep2_t _q[1];
	bn_t n;

	ep_null(_p[0]);
	ep_null(t[0]);
	ep2_null(_q[0]);
	bn_null(n);

	TRY {
		ep_new(_p[0]);
		ep_new(t[0]);
		ep2_new(_q[0]);
		bn_new(n);

		ep_norm(_p[0], p);
		ep2_norm(_q[0], q);
		ep_curve_get_ord(n);
		fp12_set_dig(r, 1);

		if (!ep_is_infty(p) && !ep2_is_infty(q)) {
			pp_mil_lit_k12(r, t, _p, _q, 1, n);
			pp_exp_k12(r, r);
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep_free(_p[0]);
		ep_free(t[0]);
		ep2_free(_q[0]);
		bn_free(n);
	}
}

void pp_map_sim_tatep_k12(fp12_t r, ep_t *p, ep2_t *q, int m) {
	ep_t *_p = RLC_ALLOCA(ep_t, m), *t = RLC_ALLOCA(ep_t, m);
	ep2_t *_q = RLC_ALLOCA(ep2_t, m);
	bn_t n;
	int i, j;

	bn_null(n);

	TRY {
		bn_new(n);
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(t[i]);
			ep2_null(_q[i]);
			ep_new(_p[i]);
			ep_new(t[i]);
			ep2_new(_q[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep2_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep2_norm(_q[j++], q[i]);
			}
		}

		ep_curve_get_ord(n);
		fp12_set_dig(r, 1);
		if (j > 0) {
			pp_mil_lit_k12(r, t, _p, _q, j, n);
			pp_exp_k12(r, r);
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(n);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep_free(t[i]);
			ep2_free(_q[i]);
		}
	}
}

#endif

#if PP_MAP == WEILP || !defined(STRIP)

void pp_map_weilp_k2(fp2_t r, ep_t p, ep_t q) {
	ep_t _p[1], _q[1], t0[1], t1[1];
	fp2_t r0, r1;
	bn_t n;

	ep_null(_p[0]);
	ep_null(_q[0]);
	ep_null(t0[0]);
	ep_null(t1[0]);
	fp2_null(r0);
	fp2_null(r1);
	bn_null(n);

	TRY {
		ep_new(_p[0]);
		ep_new(_q[0]);
		ep_new(t0[0]);
		ep_new(t1[0]);
		fp2_new(r0);
		fp2_new(r1);
		bn_new(n);

		ep_norm(_p[0], p);
		ep_norm(_q[0], q);
		ep_curve_get_ord(n);
		/* Since p has order n, we do not have to perform last iteration. */
		bn_sub_dig(n, n, 1);
		fp2_set_dig(r0, 1);
		fp2_set_dig(r1, 1);

		if (!ep_is_infty(_p[0]) && !ep_is_infty(_q[0])) {
			pp_mil_lit_k2(r0, t0, _p, _q, 1, n);
			pp_mil_k2(r1, t1, _q, _p, 1, n);
			fp2_inv(r1, r1);
			fp2_mul(r0, r0, r1);
			fp2_inv(r1, r0);
			fp2_inv_uni(r0, r0);
		}
		fp2_mul(r, r0, r1);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep_free(_p[0]);
		ep_free(_q[0]);
		ep_free(t0[0]);
		ep_free(t1[0]);
		fp2_free(r0);
		fp2_free(r1);
		bn_free(n);
	}
}

void pp_map_sim_weilp_k2(fp2_t r, ep_t *p, ep_t *q, int m) {
	ep_t *_p = RLC_ALLOCA(ep_t, m),
        *_q = RLC_ALLOCA(ep_t, m),
        *t0 = RLC_ALLOCA(ep_t, m),
        *t1 = RLC_ALLOCA(ep_t, m);
	fp2_t r0, r1;
	bn_t n;
	int i, j;

	fp2_null(r0);
	fp2_null(r1);
	bn_null(r);

	TRY {
		fp2_new(r0);
		fp2_new(r1);
		bn_new(n);
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(_q[i]);
			ep_null(t0[i]);
			ep_null(t1[i]);
			ep_new(_p[i]);
			ep_new(_q[i]);
			ep_new(t0[i]);
			ep_new(t1[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep_norm(_q[j++], q[i]);
			}
		}

		ep_curve_get_ord(n);
		bn_sub_dig(n, n, 1);
		fp2_set_dig(r0, 1);
		fp2_set_dig(r1, 1);

		if (j > 0) {
			pp_mil_lit_k2(r0, t0, _p, _q, j, n);
			pp_mil_k2(r1, t1, _q, _p, j, n);
			fp2_inv(r1, r1);
			fp2_mul(r0, r0, r1);
			fp2_inv(r1, r0);
			fp2_inv_uni(r0, r0);
		}
		fp2_mul(r, r0, r1);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp2_free(r0);
		fp2_free(r1);
		bn_free(n);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep_free(_q[i]);
			ep_free(t0[i]);
			ep_free(t1[i]);
		}
	}
}

void pp_map_weilp_k12(fp12_t r, ep_t p, ep2_t q) {
	ep_t _p[1], t0[1];
	ep2_t _q[1], t1[1];
	fp12_t r0, r1;
	bn_t n;

	ep_null(_p[0]);
	ep_null(t0[0]);
	ep2_null(_q[0]);
	ep2_null(t1[0]);
	fp12_null(r0);
	fp12_null(r1);
	bn_null(n);

	TRY {
		ep_new(_p[0]);
		ep_new(t0[0]);
		ep2_new(_q[0]);
		ep2_new(t1[0]);
		fp12_new(r0);
		fp12_new(r1);
		bn_new(n);

		ep_norm(_p[0], p);
		ep2_norm(_q[0], q);
		ep_curve_get_ord(n);
		bn_sub_dig(n, n, 1);
		fp12_set_dig(r0, 1);
		fp12_set_dig(r1, 1);

		if (!ep_is_infty(_p[0]) && !ep2_is_infty(_q[0])) {
			pp_mil_lit_k12(r0, t0, _p, _q, 1, n);
			pp_mil_k12(r1, t1, _q, _p, 1, n);
			fp12_inv(r1, r1);
			fp12_mul(r0, r0, r1);
			fp12_inv(r1, r0);
			fp12_inv_uni(r0, r0);
		}
		fp12_mul(r, r0, r1);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep_free(_p[0]);
		ep_free(t0[0]);
		ep2_free(_q[0]);
		ep2_free(t1[0]);
		fp12_free(r0);
		fp12_free(r1);
		bn_free(n);
	}
}

void pp_map_sim_weilp_k12(fp12_t r, ep_t *p, ep2_t *q, int m) {
	ep_t *_p = RLC_ALLOCA(ep_t, m), *t0 = RLC_ALLOCA(ep_t, m);
	ep2_t *_q = RLC_ALLOCA(ep2_t, m), *t1 = RLC_ALLOCA(ep2_t, m);
	fp12_t r0, r1;
	bn_t n;
	int i, j;

	fp12_null(r0);
	fp12_null(r1);
	bn_null(r);

	TRY {
		fp12_new(r0);
		fp12_new(r1);
		bn_new(n);
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(t0[i]);
			ep2_null(_q[i]);
			ep2_null(t1[i]);
			ep_new(_p[i]);
			ep_new(t0[i]);
			ep2_new(_q[i]);
			ep2_new(t1[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep2_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep2_norm(_q[j++], q[i]);
			}
		}

		ep_curve_get_ord(n);
		bn_sub_dig(n, n, 1);
		fp12_set_dig(r0, 1);
		fp12_set_dig(r1, 1);

		if (j > 0) {
			pp_mil_lit_k12(r0, t0, _p, _q, j, n);
			pp_mil_k12(r1, t1, _q, _p, j, n);
			fp12_inv(r1, r1);
			fp12_mul(r0, r0, r1);
			fp12_inv(r1, r0);
			fp12_inv_uni(r0, r0);
		}
		fp12_mul(r, r0, r1);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp12_free(r0);
		fp12_free(r1);
		bn_free(n);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep_free(t0[i]);
			ep2_free(_q[i]);
			ep2_free(t1[i]);
		}
	}
}

#endif

#if PP_MAP == OATEP || !defined(STRIP)

void pp_map_oatep_k12(fp12_t r, ep_t p, ep2_t q) {
	ep_t _p[1];
	ep2_t t[1], _q[1];
	bn_t a;

	ep_null(_p[0]);
	ep2_null(_q[0]);
	ep2_null(t[0]);
	bn_null(a);

	TRY {
		ep_new(_p[0]);
		ep2_new(_q[0]);
		ep2_new(t[0]);
		bn_new(a);

		fp_prime_get_par(a);
		fp12_set_dig(r, 1);

		ep_norm(_p[0], p);
		ep2_norm(_q[0], q);

		if (!ep_is_infty(_p[0]) && !ep2_is_infty(_q[0])) {
			switch (ep_curve_is_pairf()) {
				case EP_BN:
					bn_mul_dig(a, a, 6);
					bn_add_dig(a, a, 2);
					/* r = f_{|a|,Q}(P). */
					pp_mil_k12(r, t, _q, _p, 1, a);
					if (bn_sign(a) == RLC_NEG) {
						/* f_{-a,Q}(P) = 1/f_{a,Q}(P). */
						fp12_inv_uni(r, r);
						ep2_neg(t[0], t[0]);
					}
					pp_fin_k12_oatep(r, t[0], _q[0], _p[0]);
					pp_exp_k12(r, r);
					break;
				case EP_B12:
					/* r = f_{|a|,Q}(P). */
					pp_mil_k12(r, t, _q, _p, 1, a);
					if (bn_sign(a) == RLC_NEG) {
						fp12_inv_uni(r, r);
						ep2_neg(t[0], t[0]);
					}
					pp_exp_k12(r, r);
					break;
			}
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep_free(_p[0]);
		ep2_free(_q[0]);
		ep2_free(t[0]);
		bn_free(a);
	}
}

void pp_map_sim_oatep_k12(fp12_t r, ep_t *p, ep2_t *q, int m) {
	ep_t *_p = RLC_ALLOCA(ep_t, m);
	ep2_t* t = RLC_ALLOCA(ep2_t, m), *_q = RLC_ALLOCA(ep2_t, m);
	bn_t a;
	int i, j;

	TRY {
		bn_null(a);
		bn_new(a);
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep2_null(_q[i]);
			ep2_null(t[i]);
			ep_new(_p[i]);
			ep2_new(_q[i]);
			ep2_new(t[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep2_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep2_norm(_q[j++], q[i]);
			}
		}

		fp_prime_get_par(a);
		fp12_set_dig(r, 1);

		if (j > 0) {
			switch (ep_curve_is_pairf()) {
				case EP_BN:
					bn_mul_dig(a, a, 6);
					bn_add_dig(a, a, 2);
					/* r = f_{|a|,Q}(P). */
					pp_mil_k12(r, t, _q, _p, j, a);
					if (bn_sign(a) == RLC_NEG) {
						/* f_{-a,Q}(P) = 1/f_{a,Q}(P). */
						fp12_inv_uni(r, r);
					}
					for (i = 0; i < j; i++) {
						if (bn_sign(a) == RLC_NEG) {
							ep2_neg(t[i], t[i]);
						}
						pp_fin_k12_oatep(r, t[i], _q[i], _p[i]);
					}
					pp_exp_k12(r, r);
					break;
				case EP_B12:
					/* r = f_{|a|,Q}(P). */
					pp_mil_k12(r, t, _q, _p, j, a);
					if (bn_sign(a) == RLC_NEG) {
						fp12_inv_uni(r, r);
					}
					pp_exp_k12(r, r);
					break;
			}
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(a);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep2_free(_q[i]);
			ep2_free(t[i]);
		}
	}
}

#endif

/* degree is equal to 18 */

/**
 * Compute the Miller loop for pairings of type G_2 x G_1 over the bits of a
 * given parameter represented in sparse form.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] q				- the vector of first arguments in affine coordinates.
 * @param[in] p				- the vector of second arguments in affine coordinates.
 * @param[in] n 			- the number of pairings to evaluate.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_k18(fp18_t r, ep3_t* t, ep3_t* q, ep_t* p, int m, bn_t a) {
	fp18_t l;
	ep_t* _p = RLC_ALLOCA(ep_t, m);
	ep3_t* _q = RLC_ALLOCA(ep3_t, m);
	int i, j, len = bn_bits(a) + 1;
	int8_t* s = RLC_ALLOCA(int8_t, len);

	if (m == 0) {
		return;
	}

	fp18_null(l);

	TRY{
		fp18_new(l);
		fp18_zero(l);

		for (j = 0; j < m; j++) {
			ep_null(_p[j]);
			ep3_null(_q[j]);
			ep_new(_p[j]);
			ep3_new(_q[j]);
			ep3_copy(t[j], q[j]);
			ep3_neg(_q[j], q[j]);
			fp_add(_p[j]->x, p[j]->x, p[j]->x);
			fp_add(_p[j]->x, _p[j]->x, p[j]->x);
			fp_neg(_p[j]->y, p[j]->y);
		}
		bn_rec_naf(s, &len, a, 2);
		/*for (j = 0; j < m; j++) {
			pp_dbl_k18(l, t[j], t[j], _p[j]);
			fp18_mul_dxs(r, r, l);
		}
		if (s[len - 2] > 0) {
			for (j = 0; j < m; j++) {
				pp_add_k18(l, t[j], q[j], p[j]);
				fp18_mul_dxs(r, r, l);
			}
		}
		if (s[len - 2] < 0) {
			for (j = 0; j < m; j++) {
				pp_add_k18(l, t[j], _q[j], p[j]);
				fp18_mul_dxs(r, r, l);
			}
		}*/

		for (i = len - 2; i >= 0; i--) {
			fp18_sqr(r, r);

			for (j = 0; j < m; j++) {
				pp_dbl_k18(l, t[j], t[j], _p[j]);
				fp18_mul(r, r, l);
				if (s[i] > 0) {
					pp_add_k18(l, t[j], q[j], p[j]);
					fp18_mul(r, r, l);
				}
				if (s[i] < 0) {
					pp_add_k18(l, t[j], _q[j], p[j]);
					fp18_mul(r, r, l);
				}
			}
			/*if(i==len-3){
				fp18_print(r);
				printf("xixixixixxiix\n");
			}*/
		}
	}
		CATCH_ANY{
			THROW(ERR_CAUGHT);
	}
		FINALLY{
			fp18_free(l);
			for (j = 0; j < m; j++) {
				ep_free(_p[j]);
				ep3_free(_q[j]);
			}
	}
}


static void pp_mil_k18_3(fp18_t r, ep3_t* t, ep3_t* q, ep_t* p, int m, bn_t a) {
	fp18_t l;
	ep_t* _p = RLC_ALLOCA(ep_t, m);
	int i, j, len = bn_bits(a) + 1;
	int8_t* s = RLC_ALLOCA(int8_t, len);

	if (m == 0) {
		return;
	}

	fp18_null(l);

	TRY{
		fp18_new(l);
		fp18_zero(l);

		for (j = 0; j < m; j++) {
			ep_null(_p[j]);
			ep_new(_p[j]);
			ep3_copy(t[j], q[j]);

			fp_add(_p[j]->x, p[j]->x, p[j]->x);
			fp_add(_p[j]->x, _p[j]->x, p[j]->x);
			fp_neg(_p[j]->y, p[j]->y);
		}

		for (int i = bn_bits(a) - 2; i >= 0; i--) {
			fp18_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_k18(l, t[j], t[j], _p[j]);
				fp18_mul_dxs(r, r, l);
				if (bn_get_bit(a, i)) {
					pp_add_k18(l, t[j], q[j], p[j]);
					fp18_mul(r, r, l);
				}
			}
		}
	}
		CATCH_ANY{
			THROW(ERR_CAUGHT);
	}
		FINALLY{
			fp18_free(l);
			for (j = 0; j < m; j++) {
				ep_free(_p[j]);
			}
	}
}


/**
 * Compute the Miller loop for pairings of type G_1 x G_2 over the bits of a
 * given parameter.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] p				- the first pairing argument in affine coordinates.
 * @param[in] q				- the second pairing argument in affine coordinates.
 * @param[in] n 			- the number of pairings to evaluate.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_lit_k18(fp18_t r, ep_t* t, ep_t* p, ep3_t* q, int m, bn_t a) {
	fp18_t l;
	ep3_t* _q = RLC_ALLOCA(ep3_t, m);
	int j;

	fp18_null(l);

	TRY{
		fp18_new(l);
		for (j = 0; j < m; j++) {
			ep3_null(_q[j]);
			ep3_new(_q[j]);
			ep_copy(t[j], p[j]);
			ep3_neg(_q[j], q[j]);
		}

		fp18_zero(l);
		for (int i = bn_bits(a) - 2; i >= 0; i--) {
			fp18_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_lit_k18(l, t[j], t[j], q[j]);
				fp18_mul(r, r, l);
				if (bn_get_bit(a, i)) {
					pp_add_lit_k18(l, t[j], p[j], q[j]);
					fp18_mul(r, r, l);
				}
			}
		}
	}
		CATCH_ANY{
			THROW(ERR_CAUGHT);
	}
		FINALLY{
			fp18_free(l);
			for (j = 0; j < m; j++) {
				ep3_free(_q[j]);
			}
	}
}

/**
 * Compute the final lines for optimal ate pairings.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] q				- the first point of the pairing, in G_2.
 * @param[in] p				- the second point of the pairing, in G_1.
 * @param[in] a				- the loop parameter.
 */
static void pp_fin_k18_oatep(fp18_t r, ep3_t t, ep3_t q, ep_t p) {
	ep3_t q1;
	fp18_t tmp,l,xq1,yq1;
	bn_t three;

	fp18_null(tmp);
	fp18_null(l);
	fp18_null(xq1);
	fp18_null(yq1);
	ep3_null(q1);
	bn_null(three);

	TRY{
		ep3_new(q1);
		fp18_new(tmp);
		fp18_new(l);
		fp18_new(xq1);
		fp18_new(yq1);
		bn_new(three);

		bn_set_dig(three, 3);
		fp3_set_dig(q1->z, 1);
		fp18_set_dig(tmp,1);

		pp_mil_k18_3(tmp, q1, q, p, 1, three);
		fp18_frb(tmp, tmp, 1);
		fp18_mul(r, r, tmp);

		fp18_zero(xq1);
		fp18_zero(yq1);
		
		pp_norm_k18(q1,q1);

		ep3_frb(q1,q1,1);
		fp_copy(xq1[2][0][0], q1->x[0]);
		fp_copy(xq1[2][2][0], q1->x[1]);
		fp_copy(xq1[2][1][1], q1->x[2]);
		fp_copy(yq1[0][1][0], q1->y[0]);
		fp_copy(yq1[0][0][1], q1->y[1]);
		fp_copy(yq1[0][2][1], q1->y[2]);

		fp18_copy(l,xq1);
		
		pp_norm_k18(q1,t);

		fp_sub(l[2][0][0],q1->x[0], l[2][0][0]);
		fp_sub(l[2][2][0],q1->x[1], l[2][2][0]);
		fp_sub(l[2][1][1],q1->x[2], l[2][1][1]);
		//fp18_inv_uni(l, l);
		fp18_inv(l, l);	
		fp18_copy(tmp,yq1);
		fp_sub(tmp[0][1][0], q1->y[0], yq1[0][1][0]);
		fp_sub(tmp[0][0][1], q1->y[1], yq1[0][0][1]);
		fp_sub(tmp[0][2][1], q1->y[2], yq1[0][2][1]);
		fp18_mul_dxs(l, l, tmp);//lambda
		fp_sub(xq1[0][0][0],xq1[0][0][0],p->x);
		fp18_mul(l, l, xq1);
		fp18_neg(l,l);
		fp_sub(yq1[0][0][0], yq1[0][0][0], p->y);
		fp18_add(l, l, yq1);

		fp18_mul(r, r, l);
	} CATCH_ANY{
		THROW(ERR_CAUGHT);
	} FINALLY{
		fp18_free(tmp);
		fp18_free(l);
		fp18_free(xq1);
		fp18_free(yq1);
		ep3_free(q1);
		bn_free(three);
	}
}

void pp3_map_init(void) {
	ep3_curve_init();
}

void pp3_map_clean(void) {
	ep3_curve_clean();
}

#if PP_MAP == TATEP || !defined(STRIP)

void pp_map_tatep_k18(fp18_t r, ep_t p, ep3_t q) {
	ep_t _p[1], t[1];
	ep3_t _q[1];
	bn_t n;

	ep_null(_p[0]);
	ep_null(t[0]);
	ep3_null(_q[0]);
	bn_null(n);

	TRY{
		ep_new(_p[0]);
		ep_new(t[0]);
		ep3_new(_q[0]);
		bn_new(n);

		ep_norm(_p[0], p);
		ep3_norm(_q[0], q);
		ep_curve_get_ord(n);
		fp18_set_dig(r, 1);

		if (!ep_is_infty(p) && !ep3_is_infty(q)) {
			pp_mil_lit_k18(r, t, _p, _q, 1, n);
			pp_exp_k18(r, r);
		}
	}
		CATCH_ANY{
			THROW(ERR_CAUGHT);
	}
		FINALLY{
			ep_free(_p[0]);
			ep_free(t[0]);
			ep3_free(_q[0]);
			bn_free(n);
	}
}

void pp_map_sim_tatep_k18(fp18_t r, ep_t* p, ep3_t* q, int m) {
	ep_t* _p = RLC_ALLOCA(ep_t, m), * t = RLC_ALLOCA(ep_t, m);
	ep3_t* _q = RLC_ALLOCA(ep3_t, m);
	bn_t n;
	int i, j;

	bn_null(n);

	TRY{
		bn_new(n);
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(t[i]);
			ep3_null(_q[i]);
			ep_new(_p[i]);
			ep_new(t[i]);
			ep3_new(_q[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep3_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep3_norm(_q[j++], q[i]);
			}
		}

		ep_curve_get_ord(n);
		fp18_set_dig(r, 1);
		if (j > 0) {
			pp_mil_lit_k18(r, t, _p, _q, j, n);
			pp_exp_k18(r, r);
		}
	}
		CATCH_ANY{
			THROW(ERR_CAUGHT);
	}
		FINALLY{
			bn_free(n);
			for (i = 0; i < m; i++) {
				ep_free(_p[i]);
				ep_free(t[i]);
				ep3_free(_q[i]);
			}
	}
}

#endif

#if PP_MAP == OATEP || !defined(STRIP)

void pp_map_oatep_k18(fp18_t r, ep_t p, ep3_t q) {
	ep_t _p[1];
	ep3_t t[1], _q[1];
	bn_t a;

	ep_null(_p[0]);
	ep3_null(_q[0]);
	ep3_null(t[0]);
	bn_null(a);

	TRY{
		ep_new(_p[0]);
		ep3_new(_q[0]);
		ep3_new(t[0]);
		bn_new(a);

		fp_prime_get_par(a);
		fp18_set_dig(r, 1);

		ep_norm(_p[0], p);
		ep3_norm(_q[0], q);

		if (!ep_is_infty(_p[0]) && !ep3_is_infty(_q[0])) {
				/* r = f_{|a|,Q}(P). */
				pp_mil_k18(r, t, _q, _p, 1, a);
				if (bn_sign(a) == RLC_NEG) {
					fp18_inv_uni(r, r);
					ep3_neg(t[0], t[0]);
				}
				pp_fin_k18_oatep(r, t[0], _q[0], _p[0]);
				pp_exp_k18(r, r);

		}
	}
		CATCH_ANY{
			THROW(ERR_CAUGHT);
	}
		FINALLY{
			ep_free(_p[0]);
			ep3_free(_q[0]);
			ep3_free(t[0]);
			bn_free(a);
	}
}

void pp_map_sim_oatep_k18(fp18_t r, ep_t* p, ep3_t* q, int m) {
	ep_t* _p = RLC_ALLOCA(ep_t, m);
	ep3_t* t = RLC_ALLOCA(ep3_t, m), * _q = RLC_ALLOCA(ep3_t, m);
	bn_t a;
	int i, j;

	TRY{
		bn_null(a);
		bn_new(a);
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep3_null(_q[i]);
			ep3_null(t[i]);
			ep_new(_p[i]);
			ep3_new(_q[i]);
			ep3_new(t[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep3_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep3_norm(_q[j++], q[i]);
			}
		}

		fp_prime_get_par(a);
		fp18_set_dig(r, 1);

		if (j > 0) {
			/* r = f_{|a|,Q}(P). */
			pp_mil_k18(r, t, _q, _p, j, a);
			if (bn_sign(a) == RLC_NEG) {
				/* f_{-a,Q}(P) = 1/f_{a,Q}(P). */
				fp18_inv_uni(r, r);
			}
			for (i = 0; i < j; i++) {
				if (bn_sign(a) == RLC_NEG) {
					ep3_neg(t[i], t[i]);
				}
				pp_fin_k18_oatep(r, t[i], _q[i], _p[i]);
			}
			pp_exp_k18(r, r);
		}
	}
		CATCH_ANY{
			THROW(ERR_CAUGHT);
	}
		FINALLY{
			bn_free(a);
			for (i = 0; i < m; i++) {
				ep_free(_p[i]);
				ep3_free(_q[i]);
				ep3_free(t[i]);
			}
	}
}

#endif
