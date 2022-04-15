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
 * Tests for pairings defined over prime elliptic curves.
 *
 * @ingroup test
 */

#include <stdio.h>

#include "relic.h"
#include "relic_test.h"
#include "relic_bench.h"

unsigned long long rdtsc()
{
        unsigned long long ret;
#ifdef __x86_64__
        __asm__ volatile ("rdtsc;shlq $32,%%rdx;orq %%rdx,%%rax":"=a"(ret)::"%rdx");
#else
        __asm__ volatile ("rdtsc" : "=A" (ret));
#endif
        return ret;
}

#define xp "4B5257C16F9933101885CE46DACB939C52B6A1057B8346759871C73278EB2FFD290FA0E89CCF3\
3A4E9329D64DC5C245FDD094069D7689C46A2320249E3EA886"

#define yp "7619285FAC8D51CDC6E4C48DEDA52AE60DCA9FEB929D6284E24FC98FA746F98D21DB78AE55A32\
9AE4B5FF741169AA1731C54A75B13C4ADE1570EDAE9472CB4E"

#define xq2 "44BB61B2F53AEB4F33AAC4D34764D79ADD1CA6E9A6968A7A42B118E2889F0E20507C87996E737\
692EE045224A73549FFF5225E7558BFB7A32DB386123F0A56C"
#define xq1 "2A983B85E76A204DA4DFA41BCE87DDD4487D0C775A44DCE0051E58B22D3DFC5537EAD4E8FEF23\
7FEFCA07771D4AE338520B87450893A1C72E25BAE95001E152"
#define xq0 "6E9FF41959A5C2DEA3395BC62D796A55834B0A3E30C9753033E9D3E5AA17463C9F1D8C20E6457\
63777AC0F292D91455C88088628F9211AD86940C61CEE5D60A"
#define yq2 "3AFAEF7CC86B4512537410221A6A0808265F8D77FBBD8B8AAE1AF9721B653DD8DECB6FFEC2705\
37DA80E427499008BA3CDEB3183A4DEDDD93686740AE8ECEB5"
#define yq1 "1773104E03434A1E7E5593215AA1557CA6FD9F07AB3DEA803536FBE7902186A672C4A60383FC4\
CD5B4ED2C1476DE82B2E1ECF92A9D066D7D7B852DDAE6D573D"
#define yq0 "25E54EB6D9E8329BA8D233F665B1830784D2D640DA863DBD27A06E83A59889B4068473425E8B0\
A67CF8B19D0E9BD0B2D79BE3CB9FD76A52E8ECB279849FAE2C"
static int addition2(void) {
	int code = RLC_ERR;
	bn_t k, n;
	ep_t p, q, r, s;
	fp2_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(p);
	ep_null(q);
	ep_null(r);
	ep_null(s);
	fp2_null(e1);
	fp2_null(e2);

	TRY {
		bn_new(n);
		bn_new(k);
		ep_new(p);
		ep_new(q);
		ep_new(r);
		ep_new(s);
		fp2_new(e1);
		fp2_new(e2);

		ep_curve_get_ord(n);

		TEST_BEGIN("miller addition is correct") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			ep_copy(s, r);
			pp_add_k2(e1, r, q, p);
			pp_norm_k2(r, r);
			ep_add(s, s, q);
			ep_norm(s, s);
			TEST_ASSERT(ep_cmp(r, s) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("miller addition in affine coordinates is correct") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			ep_copy(s, r);
			fp2_zero(e1);
			fp2_zero(e2);
			pp_add_k2(e1, r, q, p);
			pp_exp_k2(e1, e1);
			pp_add_k2_basic(e2, s, q, p);
			pp_exp_k2(e2, e2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
		TEST_BEGIN("miller addition in projective coordinates is correct") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			ep_copy(s, r);
			fp2_zero(e1);
			fp2_zero(e2);
			pp_add_k2(e1, r, q, p);
			pp_exp_k2(e1, e1);
			pp_add_k2_projc(e2, s, q, p);
			pp_exp_k2(e2, e2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_BEGIN("basic projective miller addition is consistent") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			ep_copy(s, r);
			fp2_zero(e1);
			fp2_zero(e2);
			pp_add_k2_projc(e1, r, q, p);
			pp_add_k2_projc_basic(e2, s, q, p);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy-reduced projective miller addition is consistent") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			ep_copy(s, r);
			fp2_zero(e1);
			fp2_zero(e2);
			pp_add_k2_projc(e1, r, q, p);
			pp_add_k2_projc_lazyr(e2, s, q, p);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
#endif /* EP_ADD = PROJC */
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep_free(p);
	ep_free(q);
	ep_free(r);
	ep_free(s);
	fp2_free(e1);
	fp2_free(e2);
	return code;
}

static int doubling2(void) {
	int code = RLC_ERR;
	bn_t k, n;
	ep_t p, q, r, s;
	fp2_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(p);
	ep_null(q);
	ep_null(r);
	ep_null(s);
	fp2_null(e1);
	fp2_null(e2);

	TRY {
		bn_new(n);
		bn_new(k);
		ep_new(p);
		ep_new(q);
		ep_new(r);
		ep_new(s);
		fp2_new(e1);
		fp2_new(e2);

		ep_curve_get_ord(n);

		TEST_BEGIN("miller doubling is correct") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			pp_dbl_k2(e1, r, q, p);
			pp_norm_k2(r, r);
			ep_dbl(s, q);
			ep_norm(s, s);
			TEST_ASSERT(ep_cmp(r, s) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("miller doubling in affine coordinates is correct") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			fp2_zero(e1);
			fp2_zero(e2);
			pp_dbl_k2(e1, r, q, p);
			pp_exp_k2(e1, e1);
			pp_dbl_k2_basic(e2, r, q, p);
			pp_exp_k2(e2, e2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
		TEST_BEGIN("miller doubling in projective coordinates is correct") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			fp2_zero(e1);
			fp2_zero(e2);
			pp_dbl_k2(e1, r, q, p);
			pp_exp_k2(e1, e1);
			pp_dbl_k2_projc(e2, r, q, p);
			pp_exp_k2(e2, e2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_BEGIN("basic projective miller doubling is correct") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			fp2_zero(e1);
			fp2_zero(e2);
			pp_dbl_k2_projc(e1, r, q, p);
			pp_dbl_k2_projc_basic(e2, r, q, p);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy-reduced projective miller doubling is consistent") {
			ep_rand(p);
			ep_rand(q);
			ep_rand(r);
			fp2_zero(e1);
			fp2_zero(e2);
			pp_dbl_k2_projc(e1, r, q, p);
			pp_dbl_k2_projc_lazyr(e2, r, q, p);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
#endif /* EP_ADD = PROJC */
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep_free(p);
	ep_free(q);
	ep_free(r);
	ep_free(s);
	fp2_free(e1);
	fp2_free(e2);
	return code;
}

static int pairing2(void) {
	int j, code = RLC_ERR;
	bn_t k, n;
	ep_t p[2], q[2], r;
	fp2_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(r);
	fp2_null(e1);
	fp2_null(e2);

	TRY {
		bn_new(n);
		bn_new(k);
		ep_new(r);
		fp2_new(e1);
		fp2_new(e2);

		for (j = 0; j < 2; j++) {
			ep_null(p[j]);
			ep_null(q[j]);
			ep_new(p[j]);
			ep_new(q[j]);
		}

		ep_curve_get_ord(n);

		TEST_BEGIN("pairing non-degeneracy is correct") {
			ep_rand(p[0]);
			ep_rand(q[0]);
			pp_map_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep_set_infty(q[0]);
			pp_map_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("pairing is bilinear") {
			ep_rand(p[0]);
			ep_rand(q[0]);
			bn_rand_mod(k, n);
			ep_mul(r, q[0], k);
			pp_map_k2(e1, p[0], r);
			pp_map_k2(e2, p[0], q[0]);
			fp2_exp(e2, e2, k);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_k2(e2, p[0], q[0]);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_k2(e2, p[0], q[0]);
			fp2_sqr(e1, e1);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(q[0], q[0]);
			pp_map_k2(e2, p[0], q[0]);
			fp2_sqr(e1, e1);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

                TEST_BEGIN("multi-pairing is correct") {
                        ep_rand(p[i % 2]);
                        ep_rand(q[i % 2]);
                        pp_map_k2(e1, p[i % 2], q[i % 2]);
                        ep_rand(p[1 - (i % 2)]);
                        ep_set_infty(q[1 - (i % 2)]);
                        pp_map_sim_k2(e2, p, q, 2);
                        TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
                        ep_set_infty(p[1 - (i % 2)]);
                        ep_rand(q[1 - (i % 2)]);
                        pp_map_sim_k2(e2, p, q, 2);
                        TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
                        ep_set_infty(q[i % 2]);
                        pp_map_sim_k2(e2, p, q, 2);
                        TEST_ASSERT(fp2_cmp_dig(e2, 1) == RLC_EQ, end);
                        ep_rand(p[0]);
                        ep_rand(q[0]);
                        pp_map_k2(e1, p[0], q[0]);
                        ep_rand(p[1]);
                        ep_rand(q[1]);
                        pp_map_k2(e2, p[1], q[1]);
                        fp2_mul(e1, e1, e2);
                        pp_map_sim_k2(e2, p, q, 2);
                        TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
                } TEST_END;

#if PP_MAP == TATEP || PP_MAP == OATEP || !defined(STRIP)
		TEST_BEGIN("tate pairing non-degeneracy is correct") {
			ep_rand(p[0]);
			ep_rand(q[0]);
			pp_map_tatep_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_tatep_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep_set_infty(q[0]);
			pp_map_tatep_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("tate pairing is bilinear") {
			ep_rand(p[0]);
			ep_rand(q[0]);
			bn_rand_mod(k, n);
			ep_mul(r, q[0], k);
			pp_map_tatep_k2(e1, p[0], r);
			pp_map_tatep_k2(e2, p[0], q[0]);
			fp2_exp(e2, e2, k);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_tatep_k2(e2, p[0], q[0]);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_tatep_k2(e2, p[0], q[0]);
			fp2_sqr(e1, e1);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(q[0], q[0]);
			pp_map_tatep_k2(e2, p[0], q[0]);
			fp2_sqr(e1, e1);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("tate multi-pairing is correct") {
			ep_rand(p[i % 2]);
			ep_rand(q[i % 2]);
			pp_map_tatep_k2(e1, p[i % 2], q[i % 2]);
			ep_rand(p[1 - (i % 2)]);
			ep_set_infty(q[1 - (i % 2)]);
			pp_map_sim_tatep_k2(e2, p, q, 2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(p[1 - (i % 2)]);
			ep_rand(q[1 - (i % 2)]);
			pp_map_sim_tatep_k2(e2, p, q, 2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(q[i % 2]);
			pp_map_sim_tatep_k2(e2, p, q, 2);
			TEST_ASSERT(fp2_cmp_dig(e2, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep_rand(q[0]);
			pp_map_tatep_k2(e1, p[0], q[0]);
			ep_rand(p[1]);
			ep_rand(q[1]);
			pp_map_tatep_k2(e2, p[1], q[1]);
			fp2_mul(e1, e1, e2);
			pp_map_sim_tatep_k2(e2, p, q, 2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_MAP == WEIL || !defined(STRIP)
		TEST_BEGIN("weil pairing non-degeneracy is correct") {
			ep_rand(p[0]);
			ep_rand(q[0]);
			pp_map_weilp_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_weilp_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep_set_infty(q[0]);
			pp_map_weilp_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("weil pairing is bilinear") {
			ep_rand(p[0]);
			ep_rand(q[0]);
			bn_rand_mod(k, n);
			ep_mul(r, q[0], k);
			pp_map_weilp_k2(e1, p[0], r);
			pp_map_weilp_k2(e2, p[0], q[0]);
			fp2_exp(e2, e2, k);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_weilp_k2(e2, p[0], q[0]);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_weilp_k2(e2, p[0], q[0]);
			fp2_sqr(e1, e1);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(q[0], q[0]);
			pp_map_weilp_k2(e2, p[0], q[0]);
			fp2_sqr(e1, e1);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("weil multi-pairing is correct") {
			ep_rand(p[i % 2]);
			ep_rand(q[i % 2]);
			pp_map_weilp_k2(e1, p[i % 2], q[i % 2]);
			ep_rand(p[1 - (i % 2)]);
			ep_set_infty(q[1 - (i % 2)]);
			pp_map_sim_weilp_k2(e2, p, q, 2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(p[1 - (i % 2)]);
			ep_rand(q[1 - (i % 2)]);
			pp_map_sim_weilp_k2(e2, p, q, 2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(q[i % 2]);
			pp_map_sim_weilp_k2(e2, p, q, 2);
			TEST_ASSERT(fp2_cmp_dig(e2, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep_rand(q[0]);
			pp_map_weilp_k2(e1, p[0], q[0]);
			ep_rand(p[1]);
			ep_rand(q[1]);
			pp_map_weilp_k2(e2, p[1], q[1]);
			fp2_mul(e1, e1, e2);
			pp_map_sim_weilp_k2(e2, p, q, 2);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep_free(r);
	fp2_free(e1);
	fp2_free(e2);

	for (j = 0; j < 2; j++) {
		ep_free(p[j]);
		ep2_free(q[j]);
	}

        return code;
}

static int addition12(void) {
	int code = RLC_ERR;
	bn_t k, n;
	ep_t p;
	ep2_t q, r, s;
	fp12_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(p);
	ep2_null(q);
	ep2_null(r);
	ep2_null(s);
	fp12_null(e1);
	fp12_null(e2);

	TRY {
		bn_new(n);
		bn_new(k);
		ep_new(p);
		ep2_new(q);
		ep2_new(r);
		ep2_new(s);
		fp12_new(e1);
		fp12_new(e2);

		ep_curve_get_ord(n);

		TEST_BEGIN("miller addition is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			pp_add_k12(e1, r, q, p);
			pp_norm_k12(r, r);
			ep2_add(s, s, q);
			ep2_norm(s, s);
			TEST_ASSERT(ep2_cmp(r, s) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("miller addition in affine coordinates is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_add_k12(e1, r, q, p);
			pp_exp_k12(e1, e1);
			pp_add_k12_basic(e2, s, q, p);
			pp_exp_k12(e2, e2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
		TEST_BEGIN("miller addition in projective coordinates is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_add_k12(e1, r, q, p);
			pp_exp_k12(e1, e1);
			pp_add_k12_projc(e2, s, q, p);
			pp_exp_k12(e2, e2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_BEGIN("basic projective miller addition is consistent") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_add_k12_projc(e1, r, q, p);
			pp_add_k12_projc_basic(e2, s, q, p);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy-reduced projective miller addition is consistent") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_add_k12_projc(e1, r, q, p);
			pp_add_k12_projc_lazyr(e2, s, q, p);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
#endif /* EP_ADD = PROJC */
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep_free(p);
	ep2_free(q);
	ep2_free(r);
	ep2_free(s);
	fp12_free(e1);
	fp12_free(e2);
	return code;
}

static int doubling12(void) {
	int code = RLC_ERR;
	bn_t k, n;
	ep_t p;
	ep2_t q, r, s;
	fp12_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(p);
	ep2_null(q);
	ep2_null(r);
	ep2_null(s);
	fp12_null(e1);
	fp12_null(e2);

	TRY {
		bn_new(n);
		bn_new(k);
		ep_new(p);
		ep2_new(q);
		ep2_new(r);
		ep2_new(s);
		fp12_new(e1);
		fp12_new(e2);

		ep_curve_get_ord(n);

		TEST_BEGIN("miller doubling is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			pp_dbl_k12(e1, r, q, p);
			pp_norm_k12(r, r);
			ep2_dbl(s, q);
			ep2_norm(s, s);
			TEST_ASSERT(ep2_cmp(r, s) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("miller doubling in affine coordinates is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			fp12_zero(e1);
			fp12_zero(e2);
			fp_neg(p->y, p->y);
			pp_dbl_k12_basic(e2, r, q, p);
			pp_exp_k12(e2, e2);
#if EP_ADD == PROJC
			/* Precompute. */
			fp_dbl(p->z, p->x);
			fp_add(p->x, p->z, p->x);
#endif
			pp_dbl_k12(e1, r, q, p);
			pp_exp_k12(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
		TEST_BEGIN("miller doubling in projective coordinates is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			fp12_zero(e1);
			fp12_zero(e2);
			/* Precompute. */
			fp_neg(p->y, p->y);
			fp_dbl(p->z, p->x);
			fp_add(p->x, p->z, p->x);
			pp_dbl_k12_projc(e2, r, q, p);
			pp_exp_k12(e2, e2);
#if EP_ADD == BASIC
			/* Revert precomputing. */
			fp_hlv(p->x, p->z);
#endif
			pp_dbl_k12(e1, r, q, p);
			pp_exp_k12(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_BEGIN("basic projective miller doubling is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_dbl_k12_projc(e1, r, q, p);
			pp_dbl_k12_projc_basic(e2, r, q, p);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy-reduced projective miller doubling is consistent") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_dbl_k12_projc(e1, r, q, p);
			pp_dbl_k12_projc_lazyr(e2, r, q, p);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
#endif /* EP_ADD = PROJC */
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep_free(p);
	ep2_free(q);
	ep2_free(r);
	ep2_free(s);
	fp12_free(e1);
	fp12_free(e2);
	return code;
}

static void pp_mil_k18_3(fp18_t r, ep3_t* t, ep3_t* q, ep_t* p, int m, bn_t a) {
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
#if EP_ADD == BASIC
			ep_neg(_p[j], p[j]);
#else
			fp_add(_p[j]->x, p[j]->x, p[j]->x);
			fp_add(_p[j]->x, _p[j]->x, p[j]->x);
			fp_neg(_p[j]->y, p[j]->y);
#endif
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
				pp_dbl_k18_basic(l, t[j], t[j], p[j]);
				fp18_mul(r, r, l);
				if (s[i] > 0) {
					pp_add_k18_basic(l, t[j], q[j], p[j]);
					fp18_mul(r, r, l);
				}
				if (s[i] < 0) {
					pp_add_k18_basic(l, t[j], _q[j], p[j]);
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
static void pp_mil_k18_4(fp18_t r, ep3_t* t, ep3_t* q, ep_t* p, int m, bn_t a) {
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
					fp18_mul_dxs(r, r, l);
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


static int pairing12(void) {
	int j, code = RLC_ERR;
	bn_t k, n;
	ep_t p[2];
	ep2_t q[2], r;
	fp12_t e1, e2;

	bn_null(k);
	bn_null(n);
	fp12_null(e1);
	fp12_null(e2);
	ep2_null(r);

	TRY {
		bn_new(n);
		bn_new(k);
		fp12_new(e1);
		fp12_new(e2);
		ep2_new(r);

		for (j = 0; j < 2; j++) {
			ep_null(p[j]);
			ep2_null(q[j]);
			ep_new(p[j]);
			ep2_new(q[j]);
		}

		ep_curve_get_ord(n);

		TEST_BEGIN("pairing non-degeneracy is correct") {
			ep_rand(p[0]);
			ep2_rand(q[0]);
			pp_map_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep2_set_infty(q[0]);
			pp_map_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("pairing is bilinear") {
			ep_rand(p[0]);
			ep2_rand(q[0]);
			bn_rand_mod(k, n);
			ep2_mul(r, q[0], k);
			pp_map_k12(e1, p[0], r);
			pp_map_k12(e2, p[0], q[0]);
			fp12_exp(e2, e2, k);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_k12(e2, p[0], q[0]);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_k12(e2, p[0], q[0]);
			fp12_sqr(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep2_dbl(q[0], q[0]);
			pp_map_k12(e2, p[0], q[0]);
			fp12_sqr(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multi-pairing is correct") {
			ep_rand(p[i % 2]);
			ep2_rand(q[i % 2]);
			pp_map_k12(e1, p[i % 2], q[i % 2]);
			ep_rand(p[1 - (i % 2)]);
			ep2_set_infty(q[1 - (i % 2)]);
			pp_map_sim_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(p[1 - (i % 2)]);
			ep2_rand(q[1 - (i % 2)]);
			pp_map_sim_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep2_set_infty(q[i % 2]);
			pp_map_sim_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp_dig(e2, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep2_rand(q[0]);
			pp_map_k12(e1, p[0], q[0]);
			ep_rand(p[1]);
			ep2_rand(q[1]);
			pp_map_k12(e2, p[1], q[1]);
			fp12_mul(e1, e1, e2);
			pp_map_sim_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_MAP == TATEP || !defined(STRIP)
		TEST_BEGIN("tate pairing non-degeneracy is correct") {
			ep_rand(p[0]);
			ep2_rand(q[0]);
			pp_map_tatep_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_tatep_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep2_set_infty(q[0]);
			pp_map_tatep_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("tate pairing is bilinear") {
			ep_rand(p[0]);
			ep2_rand(q[0]);
			bn_rand_mod(k, n);
			ep2_mul(r, q[0], k);
			pp_map_tatep_k12(e1, p[0], r);
			pp_map_tatep_k12(e2, p[0], q[0]);
			fp12_exp(e2, e2, k);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_tatep_k12(e2, p[0], q[0]);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_tatep_k12(e2, p[0], q[0]);
			fp12_sqr(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep2_dbl(q[0], q[0]);
			pp_map_tatep_k12(e2, p[0], q[0]);
			fp12_sqr(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("tate multi-pairing is correct") {
			ep_rand(p[i % 2]);
			ep2_rand(q[i % 2]);
			pp_map_tatep_k12(e1, p[i % 2], q[i % 2]);
			ep_rand(p[1 - (i % 2)]);
			ep2_set_infty(q[1 - (i % 2)]);
			pp_map_sim_tatep_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(p[1 - (i % 2)]);
			ep2_rand(q[1 - (i % 2)]);
			pp_map_sim_tatep_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep2_set_infty(q[i % 2]);
			pp_map_sim_tatep_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp_dig(e2, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep2_rand(q[0]);
			pp_map_tatep_k12(e1, p[0], q[0]);
			ep_rand(p[1]);
			ep2_rand(q[1]);
			pp_map_tatep_k12(e2, p[1], q[1]);
			fp12_mul(e1, e1, e2);
			pp_map_sim_tatep_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_MAP == WEIL || !defined(STRIP)
		TEST_BEGIN("weil pairing non-degeneracy is correct") {
			ep_rand(p[0]);
			ep2_rand(q[0]);
			pp_map_weilp_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_weilp_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep2_set_infty(q[0]);
			pp_map_weilp_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("weil pairing is bilinear") {
			ep_rand(p[0]);
			ep2_rand(q[0]);
			bn_rand_mod(k, n);
			ep2_mul(r, q[0], k);
			pp_map_weilp_k12(e1, p[0], r);
			pp_map_weilp_k12(e2, p[0], q[0]);
			fp12_exp(e2, e2, k);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_weilp_k12(e2, p[0], q[0]);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_weilp_k12(e2, p[0], q[0]);
			fp12_sqr(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep2_dbl(q[0], q[0]);
			pp_map_weilp_k12(e2, p[0], q[0]);
			fp12_sqr(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if 0
		TEST_BEGIN("weil multi-pairing is correct") {
			ep_rand(p[i % 2]);
			ep2_rand(q[i % 2]);
			pp_map_weilp_k12(e1, p[i % 2], q[i % 2]);
			ep_rand(p[1 - (i % 2)]);
			ep2_set_infty(q[1 - (i % 2)]);
			pp_map_sim_weilp_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(p[1 - (i % 2)]);
			ep2_rand(q[1 - (i % 2)]);
			pp_map_sim_weilp_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep2_set_infty(q[i % 2]);
			pp_map_sim_weilp_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp_dig(e2, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep2_rand(q[0]);
			pp_map_weilp_k12(e1, p[0], q[0]);
			ep_rand(p[1]);
			ep2_rand(q[1]);
			pp_map_weilp_k12(e2, p[1], q[1]);
			fp12_mul(e1, e1, e2);
			pp_map_sim_weilp_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
#endif

#if PP_MAP == OATEP || !defined(STRIP)
		TEST_BEGIN("optimal ate pairing non-degeneracy is correct") {
			unsigned long long t1,t2;
			ep_rand(p[0]);
			ep2_rand(q[0]);

			TEST_ASSERT(fp12_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_oatep_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep2_set_infty(q[0]);
			pp_map_oatep_k12(e1, p[0], q[0]);
			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("optimal ate pairing is bilinear") {
			ep_rand(p[0]);
			ep2_rand(q[0]);
			bn_rand_mod(k, n);
			ep2_mul(r, q[0], k);
			pp_map_oatep_k12(e1, p[0], r);
			ep_mul(p[0], p[0], k);
			pp_map_oatep_k12(e2, p[0], q[0]);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_oatep_k12(e2, p[0], q[0]);
			fp12_sqr(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep2_dbl(q[0], q[0]);
			pp_map_oatep_k12(e2, p[0], q[0]);
			fp12_sqr(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("optimal ate multi-pairing is correct") {
			ep_rand(p[i % 2]);
			ep2_rand(q[i % 2]);
			pp_map_oatep_k12(e1, p[i % 2], q[i % 2]);
			ep_rand(p[1 - (i % 2)]);
			ep2_set_infty(q[1 - (i % 2)]);
			pp_map_sim_oatep_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(p[1 - (i % 2)]);
			ep2_rand(q[1 - (i % 2)]);
			pp_map_sim_oatep_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
			ep2_set_infty(q[i % 2]);
			pp_map_sim_oatep_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp_dig(e2, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep2_rand(q[0]);
			pp_map_oatep_k12(e1, p[0], q[0]);
			ep_rand(p[1]);
			ep2_rand(q[1]);
			pp_map_oatep_k12(e2, p[1], q[1]);
			fp12_mul(e1, e1, e2);
			pp_map_sim_oatep_k12(e2, p, q, 2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
	}
	CATCH_ANY {
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	fp12_free(e1);
	fp12_free(e2);
	ep2_free(r);

	for (j = 0; j < 2; j++) {
		ep_free(p[j]);
		ep2_free(q[j]);
	}
	return code;
}

static int addition18(void) {
	int code = RLC_ERR;
	bn_t k, n;
	ep_t p;
	ep3_t q, r, s;
	fp18_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(p);
	ep3_null(q);
	ep3_null(r);
	ep3_null(s);
	fp18_null(e1);
	fp18_null(e2);

	TRY{
		bn_new(n);
		bn_new(k);
		ep_new(p);
		ep3_new(q);
		ep3_new(r);
		ep3_new(s);
		fp18_new(e1);
		fp18_new(e2);

		ep_curve_get_ord(n);

		TEST_BEGIN("miller addition is correct") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			ep3_copy(s, r);
			pp_add_k18(e1, r, q, p);
			pp_norm_k18(r, r);
			ep3_add(s, s, q);
			ep3_norm(s, s);
			TEST_ASSERT(ep3_cmp(r, s) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("miller addition in affine coordinates is correct") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			ep3_copy(s, r);
			fp18_set_dig(e1,1);
			fp18_set_dig(e2,1);
			pp_add_k18(e1, r, q, p);
			pp_exp_k18(e1, e1);
			pp_add_k18_basic(e2, s, q, p);
			pp_exp_k18(e2, e2);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
		TEST_BEGIN("miller addition in projective coordinates is correct") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			ep3_copy(s, r);
			fp18_zero(e1);
			fp18_zero(e2);
			pp_add_k18(e1, r, q, p);
			pp_exp_k18(e1, e1);
			pp_add_k18_projc(e2, s, q, p);
			pp_exp_k18(e2, e2);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_BEGIN("basic projective miller addition is consistent") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			ep3_copy(s, r);
			fp18_zero(e1);
			fp18_zero(e2);
			pp_add_k18_projc(e1, r, q, p);
			pp_add_k18_projc_basic(e2, s, q, p);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy-reduced projective miller addition is consistent") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			ep3_copy(s, r);
			fp18_zero(e1);
			fp18_zero(e2);
			pp_add_k18_projc(e1, r, q, p);
			pp_add_k18_projc_lazyr(e2, s, q, p);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
#endif /* EP_ADD = PROJC */
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	bn_free(n);
	bn_free(k);
	ep_free(p);
	ep3_free(q);
	ep3_free(r);
	ep3_free(s);
	fp18_free(e1);
	fp18_free(e2);
	return code;
}

static int doubling18(void) {
	int code = RLC_ERR;
	bn_t k, n;
	ep_t p;
	ep3_t q, r, s;
	fp18_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(p);
	ep3_null(q);
	ep3_null(r);
	ep3_null(s);
	fp18_null(e1);
	fp18_null(e2);

	TRY{
		bn_new(n);
		bn_new(k);
		ep_new(p);
		ep3_new(q);
		ep3_new(r);
		ep3_new(s);
		fp18_new(e1);
		fp18_new(e2);

		ep_curve_get_ord(n);


		TEST_BEGIN("miller doubling is correct") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			pp_dbl_k18(e1, r, q, p);
			pp_norm_k18(r, r);
			ep3_dbl(s, q);
			ep3_norm(s, s);
			TEST_ASSERT(ep3_cmp(r, s) == RLC_EQ, end);
		} TEST_END;



#if EP_ADD == BASIC || !defined(STRIP)
		TEST_BEGIN("miller doubling in affine coordinates is correct") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			fp18_zero(e1);
			fp18_zero(e2);
			pp_dbl_k18_basic(e2, r, q, p);
			pp_exp_k18(e2, e2);
			fp_neg(p->y, p->y);
			fp_dbl(p->z, p->x);
			fp_add(p->x, p->z, p->x);
			pp_dbl_k18(e1, r, q, p);
			pp_exp_k18(e1, e1);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || !defined(STRIP)
		TEST_BEGIN("miller doubling in projective coordinates is correct") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			fp18_zero(e1);
			fp18_zero(e2);
			/* Precompute. */
			fp_neg(p->y, p->y);
			fp_dbl(p->z, p->x);
			fp_add(p->x, p->z, p->x);
			pp_dbl_k18_projc(e2, r, q, p);
			pp_exp_k18(e2, e2);
#if EP_ADD == BASIC
			/* Revert precomputing. */
			fp_hlv(p->x, p->z);
#endif
			pp_dbl_k18(e1, r, q, p);
			pp_exp_k18(e1, e1);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_BEGIN("basic projective miller doubling is correct") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			fp18_zero(e1);
			fp18_zero(e2);
			pp_dbl_k18_projc(e1, r, q, p);
			pp_dbl_k18_projc_basic(e2, r, q, p);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_BEGIN("lazy-reduced projective miller doubling is consistent") {
			ep_rand(p);
			ep3_rand(q);
			ep3_rand(r);
			fp18_zero(e1);
			fp18_zero(e2);
			pp_dbl_k18_projc(e1, r, q, p);
			pp_dbl_k18_projc_lazyr(e2, r, q, p);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
#endif /* EP_ADD = PROJC */
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	bn_free(n);
	bn_free(k);
	ep_free(p);
	ep3_free(q);
	ep3_free(r);
	ep3_free(s);
	fp18_free(e1);
	fp18_free(e2);
	return code;
}

static int pairing18(void) {
	int j, code = RLC_ERR;
	bn_t k, n;
	ep_t p[2];
	ep3_t q[2], r;
	fp18_t e1, e2;

	bn_null(k);
	bn_null(n);
	fp18_null(e1);
	fp18_null(e2);
	ep3_null(r);

	TRY{
		bn_new(n);
		bn_new(k);
		fp18_new(e1);
		fp18_new(e2);
		ep3_new(r);

		for (j = 0; j < 2; j++) {
			ep_null(p[j]);
			ep3_null(q[j]);
			ep_new(p[j]);
			ep3_new(q[j]);
		}

		ep_curve_get_ord(n);

		TEST_BEGIN("pairing non-degeneracy is correct") {
			ep_rand(p[0]);
			ep3_rand(q[0]);
			pp_map_k18(e1, p[0], q[0]);
			TEST_ASSERT(fp18_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_k18(e1, p[0], q[0]);
			TEST_ASSERT(fp18_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep3_set_infty(q[0]);
			pp_map_k18(e1, p[0], q[0]);
			TEST_ASSERT(fp18_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("pairing is bilinear") {
			ep_rand(p[0]);
			ep3_rand(q[0]);
			bn_rand_mod(k, n);
			ep3_mul_basic(r, q[0], k);
			pp_map_k18(e1, p[0], r);
			pp_map_k18(e2, p[0], q[0]);
			fp18_exp(e2, e2, k);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_k18(e2, p[0], q[0]);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_k18(e2, p[0], q[0]);
			fp18_sqr(e1, e1);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep3_dbl(q[0], q[0]);
			pp_map_k18(e2, p[0], q[0]);
			fp18_sqr(e1, e1);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("multi-pairing is correct") {
			ep_rand(p[i % 2]);
			ep3_rand(q[i % 2]);
			pp_map_k18(e1, p[i % 2], q[i % 2]);
			ep_rand(p[1 - (i % 2)]);
			ep3_set_infty(q[1 - (i % 2)]);
			pp_map_sim_k18(e2, p, q, 2);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(p[1 - (i % 2)]);
			ep3_rand(q[1 - (i % 2)]);
			pp_map_sim_k18(e2, p, q, 2);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep3_set_infty(q[i % 2]);
			pp_map_sim_k18(e2, p, q, 2);
			TEST_ASSERT(fp18_cmp_dig(e2, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep3_rand(q[0]);
			pp_map_k18(e1, p[0], q[0]);
			ep_rand(p[1]);
			ep3_rand(q[1]);
			pp_map_k18(e2, p[1], q[1]);
			fp18_mul(e1, e1, e2);
			pp_map_sim_k18(e2, p, q, 2);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_MAP == TATEP || !defined(STRIP)
		TEST_BEGIN("tate pairing non-degeneracy is correct") {
			ep_rand(p[0]);
			ep3_rand(q[0]);
			pp_map_tatep_k18(e1, p[0], q[0]);
			TEST_ASSERT(fp18_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_tatep_k18(e1, p[0], q[0]);
			TEST_ASSERT(fp18_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep3_set_infty(q[0]);
			pp_map_tatep_k18(e1, p[0], q[0]);
			TEST_ASSERT(fp18_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("tate pairing is bilinear") {
			ep_rand(p[0]);
			ep3_rand(q[0]);
			bn_rand_mod(k, n);
			ep3_mul_basic(r, q[0], k);
			pp_map_tatep_k18(e1, p[0], r);
			pp_map_tatep_k18(e2, p[0], q[0]);
			fp18_exp(e2, e2, k);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_tatep_k18(e2, p[0], q[0]);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_tatep_k18(e2, p[0], q[0]);
			fp18_sqr(e1, e1);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep3_dbl(q[0], q[0]);
			pp_map_tatep_k18(e2, p[0], q[0]);
			fp18_sqr(e1, e1);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("tate multi-pairing is correct") {
			ep_rand(p[i % 2]);
			ep3_rand(q[i % 2]);
			pp_map_tatep_k18(e1, p[i % 2], q[i % 2]);
			ep_rand(p[1 - (i % 2)]);
			ep3_set_infty(q[1 - (i % 2)]);
			pp_map_sim_tatep_k18(e2, p, q, 2);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep_set_infty(p[1 - (i % 2)]);
			ep3_rand(q[1 - (i % 2)]);
			pp_map_sim_tatep_k18(e2, p, q, 2);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
			ep3_set_infty(q[i % 2]);
			pp_map_sim_tatep_k18(e2, p, q, 2);
			TEST_ASSERT(fp18_cmp_dig(e2, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep3_rand(q[0]);
			pp_map_tatep_k18(e1, p[0], q[0]);
			ep_rand(p[1]);
			ep3_rand(q[1]);
			pp_map_tatep_k18(e2, p[1], q[1]);
			fp18_mul(e1, e1, e2);
			pp_map_sim_tatep_k18(e2, p, q, 2);
			TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

		/*
		#if PP_MAP == WEIL || !defined(STRIP)
				TEST_BEGIN("weil pairing non-degeneracy is correct") {
					ep_rand(p[0]);
					ep3_rand(q[0]);
					pp_map_weilp_k18(e1, p[0], q[0]);
					TEST_ASSERT(fp18_cmp_dig(e1, 1) != RLC_EQ, end);
					ep_set_infty(p[0]);
					pp_map_weilp_k18(e1, p[0], q[0]);
					TEST_ASSERT(fp18_cmp_dig(e1, 1) == RLC_EQ, end);
					ep_rand(p[0]);
					ep3_set_infty(q[0]);
					pp_map_weilp_k18(e1, p[0], q[0]);
					TEST_ASSERT(fp18_cmp_dig(e1, 1) == RLC_EQ, end);
				} TEST_END;

				TEST_BEGIN("weil pairing is bilinear") {
					ep_rand(p[0]);
					ep3_rand(q[0]);
					bn_rand_mod(k, n);
					ep3_mul(r, q[0], k);
					pp_map_weilp_k18(e1, p[0], r);
					pp_map_weilp_k18(e2, p[0], q[0]);
					fp18_exp(e2, e2, k);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep_mul(p[0], p[0], k);
					pp_map_weilp_k18(e2, p[0], q[0]);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep_dbl(p[0], p[0]);
					pp_map_weilp_k18(e2, p[0], q[0]);
					fp18_sqr(e1, e1);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep3_dbl(q[0], q[0]);
					pp_map_weilp_k18(e2, p[0], q[0]);
					fp18_sqr(e1, e1);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
				} TEST_END;

		#if 0
				TEST_BEGIN("weil multi-pairing is correct") {
					ep_rand(p[i % 2]);
					ep3_rand(q[i % 2]);
					pp_map_weilp_k18(e1, p[i % 2], q[i % 2]);
					ep_rand(p[1 - (i % 2)]);
					ep3_set_infty(q[1 - (i % 2)]);
					pp_map_sim_weilp_k18(e2, p, q, 2);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep_set_infty(p[1 - (i % 2)]);
					ep3_rand(q[1 - (i % 2)]);
					pp_map_sim_weilp_k18(e2, p, q, 2);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep3_set_infty(q[i % 2]);
					pp_map_sim_weilp_k18(e2, p, q, 2);
					TEST_ASSERT(fp18_cmp_dig(e2, 1) == RLC_EQ, end);
					ep_rand(p[0]);
					ep3_rand(q[0]);
					pp_map_weilp_k18(e1, p[0], q[0]);
					ep_rand(p[1]);
					ep3_rand(q[1]);
					pp_map_weilp_k18(e2, p[1], q[1]);
					fp18_mul(e1, e1, e2);
					pp_map_sim_weilp_k18(e2, p, q, 2);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
				} TEST_END;
		#endif
		#endif
		*/
		/*ep_rand(p[0]);
		ep3_rand(q[0]);
		fp_read_str(p[0]->x,xp,strlen(xp),16);
		fp_read_str(p[0]->y,yp,strlen(xp),16);
		fp_read_str(q[0]->x[0],xq0,strlen(xq0),16);
		fp_read_str(q[0]->x[1],xq1,strlen(xq1),16);
		fp_read_str(q[0]->x[2],xq2,strlen(xq2),16);
		fp_read_str(q[0]->y[0],yq0,strlen(yq0),16);
		fp_read_str(q[0]->y[1],yq1,strlen(yq1),16);
		fp_read_str(q[0]->y[2],yq2,strlen(yq2),16);
		
		bn_set_dig(k,16);
		ep3_mul_basic(r,q[0],k);
		pp_map_oatep_k18(e1, p[0], r);*/
		#if PP_MAP == OATEP || !defined(STRIP)
				TEST_BEGIN("optimal ate pairing non-degeneracy is correct") {
					unsigned long long t1,t2;
					ep_rand(p[0]);
					ep3_rand(q[0]);
					pp_map_oatep_k18(e1, p[0], q[0]);
					TEST_ASSERT(fp18_cmp_dig(e1, 1) != RLC_EQ, end);
					ep_set_infty(p[0]);
					pp_map_oatep_k18(e1, p[0], q[0]);
					TEST_ASSERT(fp18_cmp_dig(e1, 1) == RLC_EQ, end);
					ep_rand(p[0]);
					ep3_set_infty(q[0]);
					pp_map_oatep_k18(e1, p[0], q[0]);
					TEST_ASSERT(fp18_cmp_dig(e1, 1) == RLC_EQ, end);
				} TEST_END;

				TEST_BEGIN("optimal ate pairing is bilinear") {
					ep_rand(p[0]);
					ep3_rand(q[0]);
					bn_rand_mod(k, n);
					ep3_mul_basic(r, q[0], k);
					pp_map_oatep_k18(e1, p[0], r);
					ep_mul(p[0], p[0], k);
					pp_map_oatep_k18(e2, p[0], q[0]);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep_dbl(p[0], p[0]);
					pp_map_oatep_k18(e2, p[0], q[0]);
					fp18_sqr(e1, e1);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep3_dbl(q[0], q[0]);
					pp_map_oatep_k18(e2, p[0], q[0]);
					fp18_sqr(e1, e1);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
				} TEST_END;

				TEST_BEGIN("optimal ate multi-pairing is correct") {
					ep_rand(p[i % 2]);
					ep3_rand(q[i % 2]);
					pp_map_oatep_k18(e1, p[i % 2], q[i % 2]);
					ep_rand(p[1 - (i % 2)]);
					ep3_set_infty(q[1 - (i % 2)]);
					pp_map_sim_oatep_k18(e2, p, q, 2);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep_set_infty(p[1 - (i % 2)]);
					ep3_rand(q[1 - (i % 2)]);
					pp_map_sim_oatep_k18(e2, p, q, 2);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
					ep3_set_infty(q[i % 2]);
					pp_map_sim_oatep_k18(e2, p, q, 2);
					TEST_ASSERT(fp18_cmp_dig(e2, 1) == RLC_EQ, end);
					ep_rand(p[0]);
					ep3_rand(q[0]);
					pp_map_oatep_k18(e1, p[0], q[0]);
					ep_rand(p[1]);
					ep3_rand(q[1]);
					pp_map_oatep_k18(e2, p[1], q[1]);
					fp18_mul(e1, e1, e2);
					pp_map_sim_oatep_k18(e2, p, q, 2);
					TEST_ASSERT(fp18_cmp(e1, e2) == RLC_EQ, end);
				} TEST_END;
		#endif
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
end:
	bn_free(n);
	bn_free(k);
	fp18_free(e1);
	fp18_free(e2);
	ep3_free(r);

	for (j = 0; j < 2; j++) {
		ep_free(p[j]);
		ep3_free(q[j]);
	}
	return code;
}

int main(void) {
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	util_banner("Tests for the PP module", 0);

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

#if FP_PRIME==508 || FP_PRIME==676
ep3_curve_set_twist(1);
	if (ep3_curve_is_twist() == 0) {
		THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}
#else
ep2_curve_set_twist(1);
	if (ep2_curve_is_twist() == 0) {
		THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}
#endif
	ep_param_print();

	util_banner("Arithmetic", 1);

	if (ep_param_embed() == 2) {
		if (addition2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling2() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (pairing2() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	if (ep_param_embed() == 12) {
		if (addition12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling12() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (pairing12() != RLC_OK) {
			core_clean();
			return 1;
		}
	}


	if (ep_param_embed() == 18) {
		if (addition18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (doubling18() != RLC_OK) {
			core_clean();
			return 1;
		}

		if (pairing18() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	util_banner("All tests have passed.\n", 0);

	core_clean();
	return 0;
}
