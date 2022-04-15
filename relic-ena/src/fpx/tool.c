#include"ena.h"

void fp2_rdc(fp2_t c, dv2_t a) {
    fp_rdc(c[0],a[0]);
    fp_rdc(c[1],a[1]);
}

void fp3_rdc(fp3_t c, dv3_t a) {
    fp_rdc(c[0],a[0]);
    fp_rdc(c[1],a[1]);
    fp_rdc(c[2],a[2]);
}

void fp12_rdc(fp12_t c, dv12_t a) {
	for(int i=0;i<3;i++)
	{
		fp_rdc(c[0][i][0], a[0][i][0]);
		fp_rdc(c[0][i][1], a[0][i][1]);
		fp_rdc(c[1][i][0], a[1][i][0]);
		fp_rdc(c[1][i][1], a[1][i][1]);
	}
}
void fp18_rdc(fp18_t c, dv18_t a) {
	for(int i=0;i<3;i++)
	{
		fp_rdc(c[0][i][0], a[0][i][0]);
		fp_rdc(c[0][i][1], a[0][i][1]);
		fp_rdc(c[1][i][0], a[1][i][0]);
		fp_rdc(c[1][i][1], a[1][i][1]);
		fp_rdc(c[2][i][0], a[2][i][0]);
		fp_rdc(c[2][i][1], a[2][i][1]);
	}
}
void dv2_copy(dv2_t c,dv2_t a){
	dv_copy(c[0],a[0],2 * RLC_FP_DIGS);
	dv_copy(c[1],a[1],2 * RLC_FP_DIGS);
}

void fp6_muln_low(dv6_t c, fp6_t a, fp6_t b) {
	dv2_t u0, u1, u2, u3;
	fp2_t t0, t1;

	dv2_null(u0);
	dv2_null(u1);
	dv2_null(u2);
	dv2_null(u3);
	fp2_null(t0);
	fp2_null(t1);

	TRY {
		dv2_new(u0);
		dv2_new(u1);
		dv2_new(u2);
		dv2_new(u3);
		fp2_new(t0);
		fp2_new(t1);

		/* v0 = a_0b_0, v1 = a_1b_1, v2 = a_2b_2,
		 * t0 = a_1 + a_2, t1 = b_1 + b_2,
		 * u4 = u1 + u2, u5 = u0 + u1, u6 = u0 + u2 */
#ifdef RLC_FP_ROOM
		fp2_mulc_low(u0, a[0], b[0]);
		fp2_mulc_low(u1, a[1], b[1]);
		fp2_mulc_low(u2, a[2], b[2]);
		fp2_addn_low(t0, a[1], a[2]);
		fp2_addn_low(t1, b[1], b[2]);
		fp2_addd_low(c[0], u1, u2);
#else
		fp2_muln_low(u0, a[0], b[0]);
		fp2_muln_low(u1, a[1], b[1]);
		fp2_muln_low(u2, a[2], b[2]);
		fp2_addm_low(t0, a[1], a[2]);
		fp2_addm_low(t1, b[1], b[2]);
		fp2_addc_low(c[0], u1, u2);
#endif
		/* t2 (c_0) = v0 + E((a_1 + a_2)(b_1 + b_2) - v1 - v2) */
		fp2_muln_low(u3, t0, t1);
		fp2_subc_low(u3, u3, c[0]);
		fp2_norh_low(c[0], u3);
		fp2_addc_low(c[0], c[0], u0);

		/* c_1 = (a_0 + a_1)(b_0 + b_1) - v0 - v1 + Ev2 */
#ifdef RLC_FP_ROOM
		fp2_addn_low(t0, a[0], a[1]);
		fp2_addn_low(t1, b[0], b[1]);
		fp2_addd_low(c[1], u0, u1);
#else
		fp2_addm_low(t0, a[0], a[1]);
		fp2_addm_low(t1, b[0], b[1]);
		fp2_addc_low(c[1], u0, u1);
#endif
		fp2_muln_low(u3, t0, t1);
		fp2_subc_low(u3, u3, c[1]);
		fp2_norh_low(c[2], u2);
		fp2_addc_low(c[1], u3, c[2]);

		/* c_2 = (a_0 + a_2)(b_0 + b_2) - v0 + v1 - v2 */
#ifdef RLC_FP_ROOM
		fp2_addn_low(t0, a[0], a[2]);
		fp2_addn_low(t1, b[0], b[2]);
		fp2_addd_low(c[2], u0, u2);
#else
		fp2_addm_low(t0, a[0], a[2]);
		fp2_addm_low(t1, b[0], b[2]);
		fp2_addc_low(c[2], u0, u2);
#endif
		fp2_muln_low(u3, t0, t1);
		fp2_subc_low(u3, u3, c[2]);
		fp2_addc_low(c[2], u3, u1);
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		dv2_free(u0);
		dv2_free(u1);
		dv2_free(u2);
		dv2_free(u3);
		fp2_free(t0);
		fp2_free(t1);
	}
}

/*fp12 mul without modular reduction*/
void fp12_muln_low(dv12_t c, fp12_t a, fp12_t b) {
	dv6_t u0, u1, u2, u3;
	fp6_t t0, t1;

	dv6_null(u0);
	dv6_null(u1);
	dv6_null(u2);
	dv6_null(u3);
	fp6_null(t0);
	fp6_null(t1);

	TRY {
		dv6_new(u0);
		dv6_new(u1);
		dv6_new(u2);
		dv6_new(u3);
		fp6_new(t0);
		fp6_new(t1);

		/* Karatsuba algorithm. */

		/* u0 = a_0 * b_0. */
		fp6_mul_unr(u0, a[0], b[0]);
		/* u1 = a_1 * b_1. */
		fp6_mul_unr(u1, a[1], b[1]);
		/* t1 = a_0 + a_1. */
		fp6_add(t0, a[0], a[1]);
		/* t0 = b_0 + b_1. */
		fp6_add(t1, b[0], b[1]);
		/* c_1 = (a_0 + a_1) * (b_0 + b_1) */
		fp6_mul_unr(c[1], t0, t1);
		/* c_1 = c_1 - a_0b_0 - a_1b_1. */
		for (int i = 0; i < 3; i++) {
			fp2_addc_low(u3[i], u0[i], u1[i]);
			fp2_subc_low(c[1][i], c[1][i], u3[i]);
		}
		/* c_0 = a_0b_0 + v * a_1b_1. */
		fp2_nord_low(u2[0], u1[2]);
		dv_copy(u2[1][0], u1[0][0], 2 * RLC_FP_DIGS);
		dv_copy(u2[1][1], u1[0][1], 2 * RLC_FP_DIGS);
		dv_copy(u2[2][0], u1[1][0], 2 * RLC_FP_DIGS);
		dv_copy(u2[2][1], u1[1][1], 2 * RLC_FP_DIGS);
		for (int i = 0; i < 3; i++) {
			fp2_addc_low(c[0][i], u0[i], u2[i]);
		}
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		dv6_free(u0);
		dv6_free(u1);
		dv6_free(u2);
		dv6_free(u3);
		fp6_free(t0);
		fp6_free(t1);
	}
}

void fp12_sqrn_low(dv12_t c, fp12_t a) {
	fp4_t t0, t1;
	dv4_t u0, u1, u2, u3, u4;

	fp4_null(t0);
	fp4_null(t1);
	dv4_null(u0);
	dv4_null(u1);
	dv4_null(u2);
	dv4_null(u3);
	dv4_null(u4);

	TRY {
		fp4_new(t0);
		fp4_new(t1);
		dv4_new(u0);
		dv4_new(u1);
		dv4_new(u2);
		dv4_new(u3);
		dv4_new(u4);

		/* a0 = (a00, a11). */
		/* a1 = (a10, a02). */
		/* a2 = (a01, a12). */

		/* (t0,t1) = a0^2 */
		fp2_copy(t0[0], a[0][0]);
		fp2_copy(t0[1], a[1][1]);
		fp4_sqr_unr(u0, t0);

		/* (t2,t3) = 2 * a1 * a2 */
		fp2_copy(t0[0], a[1][0]);
		fp2_copy(t0[1], a[0][2]);
		fp2_copy(t1[0], a[0][1]);
		fp2_copy(t1[1], a[1][2]);
		fp4_mul_unr(u1, t0, t1);
		fp2_addc_low(u1[0], u1[0], u1[0]);
		fp2_addc_low(u1[1], u1[1], u1[1]);

		/* (t4,t5) = a2^2. */
		fp4_sqr_unr(u2, t1);

		/* c2 = a0 + a2. */
		fp2_addm_low(t1[0], a[0][0], a[0][1]);
		fp2_addm_low(t1[1], a[1][1], a[1][2]);

		/* (t6,t7) = (a0 + a2 + a1)^2. */
		fp2_addm_low(t0[0], t1[0], a[1][0]);
		fp2_addm_low(t0[1], t1[1], a[0][2]);
		fp4_sqr_unr(u3, t0);

		/* c2 = (a0 + a2 - a1)^2. */
		fp2_subm_low(t0[0], t1[0], a[1][0]);
		fp2_subm_low(t0[1], t1[1], a[0][2]);
		fp4_sqr_unr(u4, t0);

		/* c2 = (c2 + (t6,t7))/2. */
#ifdef RLC_FP_ROOM
		fp2_addd_low(u4[0], u4[0], u3[0]);
		fp2_addd_low(u4[1], u4[1], u3[1]);
#else
		fp2_addc_low(u4[0], u4[0], u3[0]);
		fp2_addc_low(u4[1], u4[1], u3[1]);
#endif
		fp_hlvd_low(u4[0][0], u4[0][0]);
		fp_hlvd_low(u4[0][1], u4[0][1]);
		fp_hlvd_low(u4[1][0], u4[1][0]);
		fp_hlvd_low(u4[1][1], u4[1][1]);

		/* (t6,t7) = (t6,t7) - c2 - (t2,t3). */
		fp2_subc_low(u3[0], u3[0], u4[0]);
		fp2_subc_low(u3[1], u3[1], u4[1]);
		fp2_subc_low(u3[0], u3[0], u1[0]);
		fp2_subc_low(u3[1], u3[1], u1[1]);

		/* c2 = c2 - (t0,t1) - (t4,t5). */
		fp2_subc_low(u4[0], u4[0], u0[0]);
		fp2_subc_low(u4[1], u4[1], u0[1]);
		fp2_subc_low(c[0][1], u4[0], u2[0]);
		fp2_subc_low(c[1][2], u4[1], u2[1]);

		/* c1 = (t6,t7) + (t4,t5) * E. */
		fp2_nord_low(u4[1], u2[1]);
		fp2_addc_low(c[1][0], u3[0], u4[1]);
		fp2_addc_low(c[0][2], u3[1], u2[0]);

		/* c0 = (t0,t1) + (t2,t3) * E. */
		fp2_nord_low(u4[1], u1[1]);
		fp2_addc_low(c[0][0], u0[0], u4[1]);
		fp2_addc_low(c[1][1], u0[1], u1[0]);

	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		fp4_free(t0);
		fp4_free(t1);
		dv4_free(u0);
		dv4_free(u1);
		dv4_free(u2);
		dv4_free(u3);
		dv4_free(u4);
	}
}

void fp18_muln_low(dv18_t c, fp18_t a, fp18_t b) {
	dv6_t u0, u1, u2, u3, u4, u5;
	fp6_t t0, t1;

	dv6_null(u0);
	dv6_null(u1);
	dv6_null(u2);
	dv6_null(u3);
	dv6_null(u4);
	dv6_null(u5);
	fp6_null(t0);
	fp6_null(t1);

	TRY {
		dv6_new(u0);
		dv6_new(u1);
		dv6_new(u2);
		dv6_new(u3);
		dv6_new(u4);
		dv6_new(u5);
		fp6_new(t0);
		fp6_new(t1);

		/* Karatsuba algorithm. */

		/* u0 = a_0 * b_0. */
		fp6_mul_unr(u0, a[0], b[0]);
		/* u1 = a_1 * b_1. */
		fp6_mul_unr(u1, a[1], b[1]);
		/* u2 = a_2 * b_2. */
		fp6_mul_unr(u2, a[2], b[2]);

		fp6_add(t0, a[1], a[2]);
		fp6_add(t1, b[1], b[2]);
		fp6_mul_unr(u3, t0, t1);
		for (int i = 0; i < 3; i++) {
			fp2_subc_low(u3[i], u3[i], u1[i]);
			fp2_subc_low(u3[i], u3[i], u2[i]);
		}
		fp2_nord_low(u4[0], u3[2]);
		fp2_addc_low(u3[2], u3[1], u0[2]);
		fp2_addc_low(u3[1], u3[0], u0[1]);
		fp2_addc_low(u3[0], u4[0], u0[0]);

		fp6_add(t0, a[0], a[1]);
		fp6_add(t1, b[0], b[1]);
		fp6_mul_unr(u4, t0, t1);
		for (int i = 0; i < 3; i++) {
			fp2_subc_low(u4[i], u4[i], u0[i]);
			fp2_subc_low(u4[i], u4[i], u1[i]);
		}
		fp2_nord_low(u5[0], u2[2]);
		dv_copy(u5[1][0], u2[0][0], 2 * RLC_FP_DIGS);
		dv_copy(u5[1][1], u2[0][1], 2 * RLC_FP_DIGS);
		dv_copy(u5[2][0], u2[1][0], 2 * RLC_FP_DIGS);
		dv_copy(u5[2][1], u2[1][1], 2 * RLC_FP_DIGS);
		for (int i = 0; i < 3; i++) {
			fp2_addc_low(u4[i], u4[i], u5[i]);
			//fp2_rdcn_low(c[1][i], u4[i]);
			dv2_copy(c[1][i], u4[i]);
		}

		fp6_add(t0, a[0], a[2]);
		fp6_add(t1, b[0], b[2]);
		fp6_mul_unr(u4, t0, t1);
		for (int i = 0; i < 3; i++) {
			fp2_subc_low(u4[i], u4[i], u0[i]);
			fp2_addc_low(u4[i], u4[i], u1[i]);
			fp2_subc_low(u4[i], u4[i], u2[i]);
			//fp2_rdcn_low(c[2][i], u4[i]);
			//fp2_rdcn_low(c[0][i], u3[i]);
			dv2_copy(c[2][i],u4[i]);
			dv2_copy(c[0][i],u3[i]);
		}
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		dv6_free(u0);
		dv6_free(u1);
		dv6_free(u2);
		dv6_free(u3);
		dv6_free(u4);
		dv6_free(u5);
		fp6_free(t0);
		fp6_free(t1);
	}
}

void fp18_sqrn_low(dv18_t c, fp18_t a) {
	dv6_t u0, u1, u2, u3, u4;
	fp6_t t0,t1;

	dv6_null(u0);
	dv6_null(u1);
	dv6_null(u2);
	dv6_null(u3);
	dv6_null(u4);
	fp6_null(t0);
	fp6_null(t1);

	TRY {
		dv6_new(u0);
		dv6_new(u1);
		dv6_new(u2);
		dv6_new(u3);
		dv6_new(u4);
		fp6_new(t0);

		/* u0 = a_0^2 */
		fp6_sqr_unr(u0, a[0]);

		/* u1 = 2 * a_1 * a_2 */
		fp6_mul_unr(u1, a[1], a[2]);
		for (int i = 0; i < 3; i++) {
			fp2_addc_low(u1[i], u1[i], u1[i]);
		}

		/* t2 = a_2^2. */
		fp6_sqr_unr(u2, a[2]);

		/* c_2 = a_0 + a_2. */
		fp6_add(t1, a[0], a[2]);

		/* u3 = (a_0 + a_2 + a_1)^2. */
		fp6_add(t0, t1, a[1]);
		fp6_sqr_unr(u3, t0);

		/* u4 = (a_0 + a_2 - a_1)^2. */
		fp6_sub(t1, t1, a[1]);
		fp6_sqr_unr(u4, t1);

		/* u4 = (c_2 + t3)/2. */
		for (int i = 0; i < 3; i++) {
			fp2_addc_low(u4[i], u4[i], u3[i]);
			fp_hlvd_low(u4[i][0], u4[i][0]);
			fp_hlvd_low(u4[i][1], u4[i][1]);
		}

		/* t3 = t3 - u4 - t1. */
		for (int i = 0; i < 3; i++) {
			fp2_subc_low(u3[i], u3[i], u4[i]);
			fp2_subc_low(u3[i], u3[i], u1[i]);
		}

		/* c_2 = u4 - t0 - t2. */
		for (int i = 0; i < 3; i++) {
			fp2_subc_low(u4[i], u4[i], u0[i]);
			fp2_subc_low(u4[i], u4[i], u2[i]);
			//fp2_rdcn_low(c[2][i], u4[i]);
			dv2_copy(c[2][i], u4[i]);
		}

		/* c_0 = t0 + t1 * E. */
		fp2_nord_low(u4[0], u1[2]);
		fp2_addc_low(u1[2], u4[0], u0[0]);
		//fp2_rdcn_low(c[0][0], u1[2]);
		dv2_copy(c[0][0], u1[2]);	
		fp2_addc_low(u1[2], u1[0], u0[1]);
		//fp2_rdcn_low(c[0][1], u1[2]);
		dv2_copy(c[0][1], u1[2]);
		fp2_addc_low(u1[2], u1[1], u0[2]);
		//fp2_rdcn_low(c[0][2], u1[2]);
		dv2_copy(c[0][2], u1[2]);

		/* c_1 = t3 + t2 * E. */
		fp2_nord_low(u4[0], u2[2]);
		fp2_addc_low(u1[0], u4[0], u3[0]);
		//fp2_rdcn_low(c[1][0], u1[0]);
		dv2_copy(c[1][0], u1[0]);
		fp2_addc_low(u1[1], u2[0], u3[1]);
		//fp2_rdcn_low(c[1][1], u1[1]);
		dv2_copy(c[1][1], u1[1]);
		fp2_addc_low(u1[2], u2[1], u3[2]);
		//fp2_rdcn_low(c[1][2], u1[2]);
		dv2_copy(c[1][2], u1[2]);
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		dv6_free(u0);
		dv6_free(u1);
		dv6_free(u2);
		dv6_free(u3);
		dv6_free(u4);
		fp6_free(t0);
		fp6_free(t1);
	}
}

void fp3_addc_low(dv3_t c,dv3_t a,dv3_t b){
	fp_addc_low(c[0],a[0],b[0]);
	fp_addc_low(c[1],a[1],b[1]);
	fp_addc_low(c[2],a[2],b[2]);
}

void fp12_addc_low(dv12_t c,dv12_t a,dv12_t b){
	for(int i=0;i<3;i++){
		fp2_addc_low(c[0][i],a[0][i],b[0][i]);
		fp2_addc_low(c[1][i],a[1][i],b[1][i]);
	}
}
void fp18_addc_low(dv18_t c,dv18_t a,dv18_t b){
	for(int i=0;i<3;i++){
		fp2_addc_low(c[0][i],a[0][i],b[0][i]);
		fp2_addc_low(c[1][i],a[1][i],b[1][i]);
		fp2_addc_low(c[2][i],a[2][i],b[2][i]);
	}
}


void fp12_subc_low(dv12_t c,dv12_t a,dv12_t b){
	for(int i=0;i<3;i++){
		fp2_subc_low(c[0][i],a[0][i],b[0][i]);
		fp2_subc_low(c[1][i],a[1][i],b[1][i]);
	}
}
void fp18_subc_low(dv18_t c,dv18_t a,dv18_t b){
	for(int i=0;i<3;i++){
		fp2_subc_low(c[0][i],a[0][i],b[0][i]);
		fp2_subc_low(c[1][i],a[1][i],b[1][i]);
		fp2_subc_low(c[2][i],a[2][i],b[2][i]);
	}
}

void fp12_dbl(fp12_t c,fp12_t a){
	fp6_dbl(c[0],a[0]);
	fp6_dbl(c[1],a[1]);
}
void fp18_dbl(fp18_t c,fp18_t a){
	fp6_dbl(c[0],a[0]);
	fp6_dbl(c[1],a[1]);
	fp6_dbl(c[2],a[2]);
}

