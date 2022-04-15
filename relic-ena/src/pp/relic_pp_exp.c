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
 * Implementation of the final exponentiation for pairings over prime curves.
 *
 * @ingroup pp
 */

#include "relic_core.h"
#include "relic_pp.h"
#include "relic_util.h"

#define exp "436E28B1125E765EEE9C990C8EA748175997EE1ACFD7DB25DFC450CB232C92873DBD1FB92AFFA\
0D35ABE53750A75678304909193DBD7D0057A09358B3D295497301158C795CA0AAD56263EB4B802\
6608BB5EBCBB42FBBCBB75C777667FC01C7A24DBAD46E9D91CDF464CC4B12A866CF7D78D1467B00\
B2DEC009E58B636909C4A4AC7C4AF4A1FBD9565E3E7FAE1B74500760A35A48DE11FDF0C06CCC7B1\
017EBFE084A26A1631F208FA8FF7B55484076D55BD3FFD51BD8E3B8A338E5751B3A10BF0F9F2572\
7D8879A5E330BA4E09C8E8A2BB76860645F25BA1AB6B928D8962C9FE06872DD7650366B41C97720\
BEEA8E520B1163DC98902B29A4A5326674621120F861C315CD7D5836EB5A045A17CB4C31629CECA\
8BDC30B86ABF6EE6DC01A7F8EFFE597E30DFE6F2FE252AA64DA2E82336EA89EE325AE06C49A627E\
0FB6DAEAEBA59995832F5DAC67167B64854A55"
/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Computes the final exponentiation of a pairing defined over a Barreto-Naehrig
 * curve.
 *
 * @param[out] c			- the result.
 * @param[in] a				- the extension field element to exponentiate.
 */
static void pp_exp_bn(fp12_t c, fp12_t a) {
	fp12_t t0, t1, t2, t3;
	bn_t x;
	const int *b;
	int l;

	fp12_null(t0);
	fp12_null(t1);
	fp12_null(t2);
	fp12_null(t3);
	bn_null(x);

	TRY {
		fp12_new(t0);
		fp12_new(t1);
		fp12_new(t2);
		fp12_new(t3);
		bn_new(x);

		/*
		 * New final exponentiation following Fuentes-Castañeda, Knapp and
		 * Rodríguez-Henríquez: Fast Hashing to G_2.
		 */
		fp_prime_get_par(x);
		b = fp_prime_get_par_sps(&l);

		/* First, compute m = f^(p^6 - 1)(p^2 + 1). */
		fp12_conv_cyc(c, a);

		/* Now compute m^((p^4 - p^2 + 1) / r). */
		/* t0 = m^2x. */
		fp12_exp_cyc_sps(t0, c, b, l);
		fp12_sqr_cyc(t0, t0);
		/* t1 = m^6x. */
		fp12_sqr_cyc(t1, t0);
		fp12_mul(t1, t1, t0);
		/* t2 = m^6x^2. */
		fp12_exp_cyc_sps(t2, t1, b, l);
		/* t3 = m^12x^3. */
		fp12_sqr_cyc(t3, t2);
		fp12_exp_cyc_sps(t3, t3, b, l);

		if (bn_sign(x) == RLC_NEG) {
			fp12_inv_uni(t0, t0);
			fp12_inv_uni(t1, t1);
			fp12_inv_uni(t3, t3);
		}

		/* t3 = a = m^12x^3 * m^6x^2 * m^6x. */
		fp12_mul(t3, t3, t2);
		fp12_mul(t3, t3, t1);

		/* t0 = b = 1/(m^2x) * t3. */
		fp12_inv_uni(t0, t0);
		fp12_mul(t0, t0, t3);

		/* Compute t2 * t3 * m * b^p * a^p^2 * [b * 1/m]^p^3. */
		fp12_mul(t2, t2, t3);
		fp12_mul(t2, t2, c);
		fp12_inv_uni(c, c);
		fp12_mul(c, c, t0);
		fp12_frb(c, c, 3);
		fp12_mul(c, c, t2);
		fp12_frb(t0, t0, 1);
		fp12_mul(c, c, t0);
		fp12_frb(t3, t3, 2);
		fp12_mul(c, c, t3);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp12_free(t0);
		fp12_free(t1);
		fp12_free(t2);
		fp12_free(t3);
		bn_free(x);
	}
}

/**
 * Computes the final exponentiation of a pairing defined over a
 * Barreto-Lynn-Scott curve.
 *
 * @param[out] c			- the result.
 * @param[in] a				- the extension field element to exponentiate.
 */
static void pp_exp_b12(fp12_t c, fp12_t a) {
	fp12_t t0, t1, t2, t3;
	bn_t x;
	const int *b;
	int l;

	fp12_null(t0);
	fp12_null(t1);
	fp12_null(t2);
	fp12_null(t3);
	bn_null(x);

	TRY {
		fp12_new(t0);
		fp12_new(t1);
		fp12_new(t2);
		fp12_new(t3);
		bn_new(x);

		/*
		 * Final exponentiation following Ghammam and Fouotsa:
		 * On the Computation of Optimal Ate Pairing at the 192-bit Level.
		 */
		fp_prime_get_par(x);
		b = fp_prime_get_par_sps(&l);

		/* First, compute m^(p^6 - 1)(p^2 + 1). */
		fp12_conv_cyc(c, a);

		/* Now compute m^((p^4 - p^2 + 1) / r). */
		/* t0 = f^2. */
		fp12_sqr_cyc(t0, c);

		/* t1 = f^x. */
		fp12_exp_cyc_sps(t1, c, b, l);
		if (bn_sign(x) == RLC_NEG) {
			fp12_inv_uni(t1, t1);
		}

		/* t2 = f^(x^2). */
		fp12_exp_cyc_sps(t2, t1, b, l);
		if (bn_sign(x) == RLC_NEG) {
			fp12_inv_uni(t2, t2);
		}

		/* t1 = t2/(t1^2 * f). */
		fp12_inv_uni(t3, c);
		fp12_sqr_cyc(t1, t1);
		fp12_mul(t1, t1, t3);
		fp12_inv_uni(t1, t1);
		fp12_mul(t1, t1, t2);

		/* t2 = t1^x. */
		fp12_exp_cyc_sps(t2, t1, b, l);
		if (bn_sign(x) == RLC_NEG) {
			fp12_inv_uni(t2, t2);
		}

		/* t3 = t2^x/t1. */
		fp12_exp_cyc_sps(t3, t2, b, l);
		if (bn_sign(x) == RLC_NEG) {
			fp12_inv_uni(t3, t3);
		}
		fp12_inv_uni(t1, t1);
		fp12_mul(t3, t1, t3);

		/* t1 = t1^(-p^3 ) * t2^(p^2). */
		fp12_inv_uni(t1, t1);
		fp12_frb(t1, t1, 3);
		fp12_frb(t2, t2, 2);
		fp12_mul(t1, t1, t2);

		/* t2 = f * f^2 * t3^x. */
		fp12_exp_cyc_sps(t2, t3, b, l);
		if (bn_sign(x) == RLC_NEG) {
			fp12_inv_uni(t2, t2);
		}
		fp12_mul(t2, t2, t0);
		fp12_mul(t2, t2, c);

		/* Compute t1 * t2 * t3^p. */
		fp12_mul(t1, t1, t2);
		fp12_frb(t2, t3, 1);
		fp12_mul(c, t1, t2);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp12_free(t0);
		fp12_free(t1);
		fp12_free(t2);
		fp12_free(t3);
		bn_free(x);
	}
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void pp_exp_k2(fp2_t c, fp2_t a) {
	bn_t e, n;

	bn_null(n);
	bn_null(e);

	TRY {
		bn_new(n);
		bn_new(e);

		ep_curve_get_ord(n);

		fp2_conv_uni(c, a);
		dv_copy(e->dp, fp_prime_get(), RLC_FP_DIGS);
		e->used = RLC_FP_DIGS;
		e->sign = RLC_POS;
		bn_add_dig(e, e, 1);
		bn_div(e, e, n);
		fp2_exp_uni(c, c, e);
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		bn_free(n);
		bn_free(e);
	}
}

void pp_exp_k12(fp12_t c, fp12_t a) {
	switch (ep_curve_is_pairf()) {
		case EP_BN:
			pp_exp_bn(c, a);
			break;
		case EP_B12:
			pp_exp_b12(c, a);
			break;
	}
}

void pp_exp_k18(fp18_t c, fp18_t a) {
	fp18_t t[10],t0,t1,t2,t3,t4,t5,t6;
	bn_t x;
	const int* b;
	int l;

	for (int i = 0; i < 10; i++)fp18_null(t[i]);
	fp18_null(t0);
	fp18_null(t1);
	fp18_null(t2);
	fp18_null(t3);
	fp18_null(t4);
	fp18_null(t5);
	fp18_null(t6);
	bn_null(x);

	TRY{
		for (int i = 0; i < 10; i++)fp18_new(t[i]);
		fp18_new(t0);
		fp18_new(t1);
		fp18_new(t2);
		fp18_new(t3);
		fp18_new(t4);
		fp18_new(t5);
		fp18_new(t6);
		bn_new(x);

		fp_prime_get_par(x);
		b = fp_prime_get_par_sps(&l);

		//First, compute m = f^(p^9 - 1)(p^3 + 1).
		fp18_conv_cyc(a, a);

		fp18_inv(t0, a);//a=1/f
		fp18_frb(t[0], t0, 1);//y23
		fp18_sqr_cyc(t[0], t[0]);
		fp18_frb(t[1], t0, 4);//y21
		fp18_mul(t[0], t[0], t[1]);
		fp18_exp_cyc_sps(t0, a, b, l);//t0=f^-x
		fp18_inv_uni(t1, t0);//t1=f^x
		fp18_frb(t[8], t0, 1);
		fp18_mul(t[0], t[0], t[8]);
		fp18_frb(t[1], t1, 3);
		fp18_mul(t[1], t[1], t[0]);
		fp18_exp_cyc_sps(t2, t1, b, l);//t2=f^-x^2
		fp18_frb(t[2], t2, 1);
		fp18_mul(t[0], t[2], t[0]);
		fp18_exp_cyc_sps(t3, t2, b, l);//t3=f^x^3
		fp18_inv_uni(t3, t3);//t3=f^-x^3
		fp18_frb(t[2], t3, 2);
		fp18_mul(t[0], t[2], t[0]);
		fp18_frb(t[3], t0, 4);
		fp18_frb(t[6], a, 2);
		fp18_mul(t[6],t[6],t[3]);
		fp18_mul(t[1], t[1], t[3]);
		fp18_exp_cyc_sps(t5, t3, b, l);//t5=f^x^4
		fp18_exp_cyc_sps(t4, t5, b, l);//t4=1/f^x^5
		fp18_inv_uni(t0, t5);
		fp18_frb(t[2], t4, 4);
		fp18_frb(t[4], a, 5);
		fp18_mul(t[2], t[2], t[4]);//y2
		fp18_mul(t[7], t[2], t[3]);
		fp18_mul(t[2], t[1], t1);
		fp18_mul(t[0], t[0], t[1]);//y17
		fp18_inv_uni(t1, t2);//t1=f^x^2
		fp18_frb(t[3],t1, 3);//y18
		fp18_mul(t[1], t[3], t[8]);
		fp18_frb(t[5], t0, 2);//y9
		fp18_mul(t[3], t[3], t[5]);
		fp18_mul(t[2], t[1], t[2]);
		fp18_frb(t[4], t0, 4);//y8
		fp18_mul(t[1], t[4], t[1]);
		fp18_frb(t[5], t2, 2);//y16
		fp18_mul(t[5], t[5], t[2]);
		fp18_frb(t[8], t0, 1);
		fp18_frb(t6, t5, 3);//y7
		fp18_mul(t[2], t[2], t6);
		fp18_mul(t[4], t[5], t1);//y15
		fp18_mul(t[5], t[5], t[8]);//y11
		fp18_mul(t[0], t[4], t[0]);
		fp18_inv_uni(t6, t4);
		fp18_frb(t0, t6, 3);//y6
		fp18_mul(t[4], t[4], t0);
		fp18_sqr_cyc(t[0], t[0]);
		fp18_frb(t0, t2, 4);//y13
		fp18_mul(t[0], t[0], t0);
		fp18_frb(t[8], t3, 1);
		fp18_mul(t[8], t[0], t[8]);
		fp18_inv_uni(t0, t3);
		fp18_frb(t[9], t0, 3);
		fp18_mul(t[8], t[8], t[9]);//y12
		fp18_mul(t[0], t[2], t[0]);
		fp18_mul(t[2], t[8], t[5]);
		fp18_mul(t[5], t[4], t[8]);
		fp18_mul(t[1], t[2], t[1]);
		fp18_mul(t[2], t[2], t[5]);
		fp18_frb(t[9], t3, 4);
		fp18_mul(t[8], t[9], t0);//y10
		fp18_mul(t[3], t[3], t[8]);
		fp18_exp_cyc_sps(t3, t4, b, l);//t3=f^x^6
		fp18_frb(t[4], t3, 3);//y1
		fp18_mul(t[4], t[8], t[4]);
		fp18_mul(t[0], t[0], t[3]);
		fp18_frb(t[9], t0, 5);
		fp18_mul(t[9], t[9], t6);
		fp18_frb(t[8], t4, 1);
		fp18_mul(t[8], t[9], t[8]);
		fp18_inv_uni(t6, t3);
		fp18_frb(t[9], t6, 2);
		fp18_mul(t[8], t[9], t[8]);
		fp18_mul(t[3], t[3], t[8]);
		fp18_mul(t[0], t[1], t[0]);//y3
		fp18_frb(t[9], t5, 5);
		fp18_mul(t[8], t[9], t3);
		fp18_exp_cyc_sps(t3, t3, b, l);//t3=1/f^x^7
		fp18_frb(t[9], t3, 2);
		fp18_mul(t[1], t[9], t[1]);
		fp18_mul(t[1], t[1], t[8]);//y0
		fp18_sqr_cyc(t[0], t[0]);
		fp18_frb(t[9], t4, 2);
		fp18_mul(t[0], t[0], t5);
		fp18_mul(t[0], t[0], t[9]);
		fp18_frb(t[9], t1, 5);
		fp18_mul(t[0], t[0], t[9]);//y5
		fp18_sqr_cyc(t[2], t[2]);
		fp18_mul(t[2], t[2], t[3]);
		fp18_mul(t[3], t[0], t[6]);
		fp18_mul(t[0], t[0], t[1]);
		fp18_sqr_cyc(t[1], t[3]);
		fp18_mul(t[1], t[1], t[2]);
		fp18_mul(t[1], t[1], t[7]);
		fp18_mul(t[2], t[2], t[4]);
		fp18_mul(t[2], t[2], t[1]);
		fp18_mul(t[0], t[0], t[1]);
		fp18_sqr_cyc(c, t[2]);
		fp18_mul(c, c, t[0]);
	}
		CATCH_ANY{
			THROW(ERR_CAUGHT);
	}
		FINALLY{
			for (int i = 0; i < 10; i++)fp18_free(t[i]);
			fp18_free(t0);
			fp18_free(t1);
			fp18_free(t2);
			fp18_free(t3);
			fp18_free(t4);
			fp18_free(t5);
			fp18_free(t6);
			bn_free(x);
	}

}

