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
 * Implementation of configuration of prime elliptic curves over quadratic
 * extensions.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if defined(EP_ENDOM) && FP_PRIME == 508
/**
 * Parameters for a KSS18 curve over a cubic extension.
 */
/** @{ */
#define KSS_P508_A0		"0"
#define KSS_P508_A1		"0"
#define KSS_P508_A2		"0"
#define KSS_P508_B0		"0"
#define KSS_P508_B1		"0"
#define KSS_P508_B2		"C33B72B87E5A44DF420B8A16E4727D4A5524C027B6231A491E008FC121F07073D3F6EB93785360E0B1B3D5304D0957E6B5CC3A8F69C137ACD1CCCF628BC1B1C"
#define KSS_P508_X0		"3586D32255A63EA229C176C6760C6DFC5DC81326AF6D61EC70B49BF4A9C43CB8EDD3EACED364A127113F43BB12EACE094016A02BBB5A1E718D3516F6470F510"
#define KSS_P508_X1		"8A4BD501BB0D6C0C9DB8ADE27CC596A667D92EC26028FA38CE4972A0C02EF6C1E60D95FE8A214352FDAB41F72EFF81ED8164A18393F221BB575286300342D2E"
#define KSS_P508_X2		"9D46B10306EDB4C77FA4254F292D4E85B7A8141B8B46D1C777947E109BFFA9E48977B7C7C26BC331E81B75275E470E3FB2C81ECA3224597491E5AD9E32C3F3E"
#define KSS_P508_Y0		"5BC1AEC8ABB89B5C1CEE28F6FBE63E3CAD4039C040CC195F3ED0747941AFA45727984158D0D15C1DB912BCDA3C1ABAB388B9D3B055E208A4460C0B3D84CCF77"
#define KSS_P508_Y1		"797C889D8B43DCAC840B7DB6ADEF6D1798FD2FCC101503A807F8B8B80446C947BF171C59CEC13CFF4E73D6339EB88CE8800BC4F294570F28DF5E3E328742B4"
#define KSS_P508_Y2		"1BBC8D443860C0A2DAA2CE2EAF980FA576F9CA1E92F5CB831CC6C4312563DA859E04A715DDB7CDF334BD0402D309487B0146539263E311A2E2AD68A2E9A2D43"
#define KSS_P508_R		"BF33E1C9934E7868ECE51D291E5644DA8A2F179CEE74854EE6819B240F20CE4E7D19F4CDABA6EAEA5B0E3000000001"
/** @} */
#endif

#if defined(EP_ENDOM) && FP_PRIME ==676
/**
 * Parameters for a KSS18 curve over a cubic extension.
 */
/** @{ */
#define KSS_P676_A0		"0"
#define KSS_P676_A1		"0"
#define KSS_P676_A2		"0"
#define KSS_P676_B0		"0"
#define KSS_P676_B1		"0"
#define KSS_P676_B2		"C30C30C30C30DC30C2FFE18619CC8616FE1618C43A2DA5ACE3BF81206305E834BC180EE15B76C\
D748BD5DCFBF65CCC35F527C35E771FBD43616D4719941AAB02643C88002A0A2030CF17D0E587723453FF33EA0DC"
#define KSS_P676_X0		"88CE5C16E857C77841753C1AA9E2B84FBD9CF8C48950C21059503E85F5625828AAC4F9A9DBE40\
E39E3DE3A76698490263B84A1B9C4822F5BB9FB6E5918CBBE72D52D2B56001E7D924EB3113A01A1DA2D82D607983"
#define KSS_P676_X1		"6326813F8A29E1A73FFA110FF1DF685D40EC26AA4592F9E8D13AC881532FC4E518C41B6E46A1E\
107E0FD19E152F31F816DCA1A97F395081CAB190C58BAADF75F2C58EF8BB73EF79CD3FBABB019C9AECC614A2E57C"
#define KSS_P676_X2		"A2AF998B3E6E3713F6997D7A777AD0E1B39385886624E5C52C61921A2100917F8AEDD461E8554\
D8F5864FD10AE10CB96353F72DF3A9568C49DDBE2BE5B790057FB20E76A9FD6EC9317B107667215C3F27CC78228F"
#define KSS_P676_Y0		"4573A6FD97545C068B91A4A7425A48BDA933B2CA8F9D3E8F8BFA8C45369664A2ED95C1FF6DECF\
6A326ACAEA8DBFE533D8F8BC1B762C16934E6327A4734F032A7A26568D5ABFF18BDE9D38DBAF3B2E65DE8371FEED"
#define KSS_P676_Y1		"426A754E7790F66906263677E00F3F1DBC0807B3D420B0935DE942317977269B74A88C258407A\
F22F7C5C671C5C48068AB3A433D6C349920256297A8D313F5065AB6644389D0F2ADF9DBEDAFD74367BA03EEC2CF0"
#define KSS_P676_Y2		"5869FB3C6597E81FC6A2BC7F4451E488B782E267FE17DE54915D7EF6B165E5980E0D700E28CE3\
0E95FFB6B13FA206B0D1E5C1BEBAE014E32E9EB8153536B7DDF5ED90DD52DFE36BC56693B5601713B313AD42EB92"
#define KSS_P676_R		"2FC44AA2B49E3ED5752B49E3A3A34F1CED6F3509FABC72F1782A87A4C0981692AD2959E9865E24242D3CE8B045AF49DF91E16BDFDD482F4CD86EC6C38C0001"
/** @} */
#endif

/**
 * Assigns a set of ordinary elliptic curve parameters.
 *
 * @param[in] CURVE		- the curve parameters to assign.
 */
#define ASSIGN(CURVE)														\
	RLC_GET(str, CURVE##_A0, sizeof(CURVE##_A0));							\
	fp_read_str(a[0], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A1, sizeof(CURVE##_A1));							\
	fp_read_str(a[1], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A2, sizeof(CURVE##_A2));							\
	fp_read_str(a[2], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B0, sizeof(CURVE##_B0));							\
	fp_read_str(b[0], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B1, sizeof(CURVE##_B1));							\
	fp_read_str(b[1], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B2, sizeof(CURVE##_B2));							\
	fp_read_str(b[2], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X0, sizeof(CURVE##_X0));							\
	fp_read_str(g->x[0], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X1, sizeof(CURVE##_X1));							\
	fp_read_str(g->x[1], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X2, sizeof(CURVE##_X2));							\
	fp_read_str(g->x[2], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y0, sizeof(CURVE##_Y0));							\
	fp_read_str(g->y[0], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y1, sizeof(CURVE##_Y1));							\
	fp_read_str(g->y[1], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y2, sizeof(CURVE##_Y2));							\
	fp_read_str(g->y[2], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_R, sizeof(CURVE##_R));								\
	bn_read_str(r, str, strlen(str), 16);									\

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void ep3_curve_init(void) {
	ctx_t *ctx = core_get();

#ifdef EP_PRECO
	for (int i = 0; i < RLC_EP_TABLE; i++) {
		ctx->ep3_ptr[i] = &(ctx->ep3_pre[i]);
	}
#endif

#if ALLOC == DYNAMIC || ALLOC == STACK
	ctx->ep3_g.x[0] = ctx->ep3_gx[0];
	ctx->ep3_g.x[1] = ctx->ep3_gx[1];
	ctx->ep3_g.y[0] = ctx->ep3_gy[0];
	ctx->ep3_g.y[1] = ctx->ep3_gy[1];
	ctx->ep3_g.z[0] = ctx->ep3_gz[0];
	ctx->ep3_g.z[1] = ctx->ep3_gz[1];
#endif

#ifdef EP_PRECO
#if ALLOC == DYNAMIC
	for (int i = 0; i < RLC_EP_TABLE; i++) {
		fp3_new(ctx->ep3_pre[i].x);
		fp3_new(ctx->ep3_pre[i].y);
		fp3_new(ctx->ep3_pre[i].z);
	}
#elif ALLOC == STACK
	for (int i = 0; i < RLC_EP_TABLE; i++) {
		ctx->ep3_pre[i].x[0] = ctx->_ep3_pre[3 * i][0];
		ctx->ep3_pre[i].x[1] = ctx->_ep3_pre[3 * i][1];
		ctx->ep3_pre[i].y[0] = ctx->_ep3_pre[3 * i + 1][0];
		ctx->ep3_pre[i].y[1] = ctx->_ep3_pre[3 * i + 1][1];
		ctx->ep3_pre[i].z[0] = ctx->_ep3_pre[3 * i + 2][0];
		ctx->ep3_pre[i].z[1] = ctx->_ep3_pre[3 * i + 2][1];
	}
#endif
#endif
	ep3_set_infty(&(ctx->ep3_g));
	bn_init(&(ctx->ep3_r), RLC_FP_DIGS);
	bn_init(&(ctx->ep3_h), RLC_FP_DIGS);
}

void ep3_curve_clean(void) {
	ctx_t *ctx = core_get();
#ifdef EP_PRECO
	for (int i = 0; i < RLC_EP_TABLE; i++) {
		fp3_free(ctx->ep3_pre[i].x);
		fp3_free(ctx->ep3_pre[i].y);
		fp3_free(ctx->ep3_pre[i].z);
	}
#endif
	bn_clean(&(ctx->ep3_r));
	bn_clean(&(ctx->ep3_h));
}

int ep3_curve_is_twist(void) {
	return core_get()->ep3_is_twist;
}

void ep3_curve_get_gen(ep3_t g) {
	ep3_copy(g, &(core_get()->ep3_g));
}

void ep3_curve_get_a(fp3_t a) {
	ctx_t *ctx = core_get();
	fp_copy(a[0], ctx->ep3_a[0]);
	fp_copy(a[1], ctx->ep3_a[1]);
	fp_copy(a[2], ctx->ep3_a[2]);
}

void ep3_curve_get_b(fp3_t b) {
	ctx_t *ctx = core_get();
	fp_copy(b[0], ctx->ep3_b[0]);
	fp_copy(b[1], ctx->ep3_b[1]);
	fp_copy(b[2], ctx->ep3_b[2]);
}
//???????
void ep3_curve_get_vs(bn_t *v) {
	bn_t x, t;

	bn_null(x);
	bn_null(t);

	TRY {
		bn_new(x);
		bn_new(t);

		fp_prime_get_par(x);
		bn_copy(v[1], x);
		bn_copy(v[2], x);
		bn_copy(v[3], x);

		/* t = 2x^2. */
		bn_sqr(t, x);
		bn_dbl(t, t);

		/* v0 = 2x^2 + 3x + 1. */
		bn_mul_dig(v[0], x, 3);
		bn_add_dig(v[0], v[0], 1);
		bn_add(v[0], v[0], t);

		/* v3 = -(2x^2 + x). */
		bn_add(v[3], v[3], t);
		bn_neg(v[3], v[3]);

		/* v1 = 12x^3 + 8x^2 + x, v2 = 6x^3 + 4x^2 + x. */
		bn_dbl(t, t);
		bn_add(v[2], v[2], t);
		bn_dbl(t, t);
		bn_add(v[1], v[1], t);
		bn_rsh(t, t, 2);
		bn_mul(t, t, x);
		bn_mul_dig(t, t, 3);
		bn_add(v[2], v[2], t);
		bn_dbl(t, t);
		bn_add(v[1], v[1], t);
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		bn_free(x);
		bn_free(t);
	}
}

void ep3_curve_get_ord(bn_t n) {
	ctx_t *ctx = core_get();
	if (ctx->ep3_is_twist) {
		ep_curve_get_ord(n);
	} else {
		bn_copy(n, &(ctx->ep3_r));
	}
}

void ep3_curve_get_cof(bn_t h) {
	bn_copy(h, &(core_get()->ep3_h));
}

#if defined(EP_PRECO)

ep3_t *ep3_curve_get_tab(void) {
#if ALLOC == AUTO
	return (ep3_t *)*(core_get()->ep3_ptr);
#else
	return core_get()->ep3_ptr;
#endif
}

#endif

void ep3_curve_set_twist(int type) {
	char str[4 * RLC_FP_BYTES + 1];
	ctx_t *ctx = core_get();
	ep3_t g;
	fp3_t a;
	fp3_t b;
	bn_t r;

	ep3_null(g);
	fp3_null(a);
	fp3_null(b);
	bn_null(r);

	ctx->ep3_is_twist = 0;
	if (type == EP_MTYPE || type == EP_DTYPE) {
		ctx->ep3_is_twist = type;
	} else {
		return;
	}

	TRY {
		ep3_new(g);
		fp3_new(a);
		fp3_new(b);
		bn_new(r);

		switch (ep_param_get()) {
#if FP_PRIME ==508
			//printf("??????????");
			case KSS_P508:
				ASSIGN(KSS_P508);
				break;
#elif FP_PRIME ==676
			case KSS_P676:
				ASSIGN(KSS_P676);
				break;
#endif
			default:
				(void)str;
				THROW(ERR_NO_VALID);
				break;
		}
		fp3_zero(g->z);
		fp_set_dig(g->z[0], 1);
		g->norm = 1;

		ep3_copy(&(ctx->ep3_g), g);
		fp_copy(ctx->ep3_a[0], a[0]);
		fp_copy(ctx->ep3_a[1], a[1]);
		fp_copy(ctx->ep3_a[2], a[2]);
		fp_copy(ctx->ep3_b[0], b[0]);
		fp_copy(ctx->ep3_b[1], b[1]);
		fp_copy(ctx->ep3_b[2], b[2]);
		bn_copy(&(ctx->ep3_r), r);
		/* TODO: support non-trivial cofactors. */
		bn_set_dig(&(ctx->ep3_h), 1);
		/* I don't have a better place for this. */
		fp_prime_calc();

/*#if defined(EP_PRECO)
		ep3_mul_pre((ep3_t *)ep3_curve_get_tab(), &(ctx->ep3_g));
#endif*/
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep3_free(g);
		fp3_free(a);
		fp3_free(b);
		bn_free(r);
	}
}

void ep3_curve_set(fp3_t a, fp3_t b, ep3_t g, bn_t r, bn_t h) {
	ctx_t *ctx = core_get();
	ctx->ep3_is_twist = 0;

	fp_copy(ctx->ep3_a[0], a[0]);
	fp_copy(ctx->ep3_a[1], a[1]);
	fp_copy(ctx->ep3_a[2], a[2]);
	fp_copy(ctx->ep3_b[0], b[0]);
	fp_copy(ctx->ep3_b[1], b[1]);
	fp_copy(ctx->ep3_b[2], b[2]);

	ep3_norm(&(ctx->ep3_g), g);
	bn_copy(&(ctx->ep3_r), r);
	bn_copy(&(ctx->ep3_h), h);

#if defined(EP_PRECO)
	ep3_mul_pre((ep3_t *)ep3_curve_get_tab(), &(ctx->ep3_g));
#endif
}
