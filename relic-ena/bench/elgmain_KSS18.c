#include"ena.h"
#include"time.h"
//rdc rule:for add and sub, the two element should be >p or <p both;for mul and sqr,the two element should be <p both
/*KSS18-P676 curve parameter:
//y^2 = x^3 + 2
A=0;
B=2;
bata=-2;
Fq^2=Fq[u]/(u^2-beta);
ksai=u;
Fq^6=Fq^2[v]/(v^3-ksai);
Fq^18=Fp^6[w]/(w^3-v);
*/
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

/**
 * @file
 *
 * Benchmarks for the optimal ate pairing defined over the KSS18-P676 curve (elliptic net algorithm).
 *
 * @ingroup bench
 */

static int pairing(){
	int code = RLC_ERR;
    	fp_t B1;
	fp3_t A,B;
	fp18_t F1,F2;
	bn_t r,s,x;
	ep_t P;
	ep3_t Q,R;

  	fp_null(B1);
	fp3_null(A);
	fp3_null(B);
	bn_null(x);
	bn_null(r);
	bn_null(s);
	fp18_null(F1);
	fp18_null(F2);
	ep_null(P);
	ep3_null(Q);
	ep3_null(R);

	TRY{
	fp_new(B1);
	fp3_new(A);
	fp3_new(B);
	bn_new(x);
	bn_new(r);
	bn_new(s);
	fp18_new(F1);
	fp18_new(F2);
	ep_new(P);
	ep3_new(Q);
	ep3_new(R);

	ep_curve_get_ord(r);
	fp_prime_get_par(x);
	ep3_curve_get_a(A);//A=0
	ep3_curve_get_b(B);//B=-u^2
	fp_set_dig(B1,2);

/*==========================================================================*/
/*                                  Testing                                 */
/*==========================================================================*/
	printf("The correctness of the elliptic net algorithm on the untwist KSS curve:\n");

	printf("\nTesting if optimal ate pairing non-degeneracy is correct:\n\n");

	TEST_BEGIN("The Elliptic Net Algorithm on KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_kss(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_kss(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_kss(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The Elliptic Net Algorithm with Lazy Reduction on KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_kss_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_kss_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_kss_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm on KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_cz_kss(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_kss(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_cz_kss(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_cz_kss_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_kss_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_cz_kss_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv) on KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_cz_ri_kss(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_ri_kss(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_cz_ri_kss(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction on KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_cz_ri_kss_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_ri_kss_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_cz_ri_kss_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	printf("\nTesting if optimal ate pairing is bilinear:\n\n");

	TEST_BEGIN("The Elliptic Net Algorithm on KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_kss(F1,R,P,A[0],B1,x,1);
		opt_en_kss(F2,Q,P,A[0],B1,x,1);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_kss(F2,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The Elliptic Net Algorithm with Lazy Reduction on KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_kss_lz(F1,R,P,A[0],B1,x,2);
		opt_en_kss_lz(F2,Q,P,A[0],B1,x,2);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_kss_lz(F2,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm on KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_cz_kss(F1,R,P,A[0],B1,x,3);
		opt_en_cz_kss(F2,Q,P,A[0],B1,x,3);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_kss(F2,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_cz_kss_lz(F1,R,P,A[0],B1,x,4);
		opt_en_cz_kss_lz(F2,Q,P,A[0],B1,x,4);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_kss_lz(F2,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv) on KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_cz_ri_kss(F1,R,P,A[0],B1,x,1);
		opt_en_cz_ri_kss(F2,Q,P,A[0],B1,x,1);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_ri_kss(F2,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction on KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_cz_ri_kss_lz(F1,R,P,A[0],B1,x,2);
		opt_en_cz_ri_kss_lz(F2,Q,P,A[0],B1,x,2);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_ri_kss_lz(F2,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	printf("\nThe correctness of the elliptic net algorithm on the twisted KSS curve:\n");

	printf("\nTesting if optimal ate pairing non-degeneracy is correct:\n\n");

	TEST_BEGIN("The Elliptic Net Algorithm on twsited KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_tw_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_tw_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_tw_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The Elliptic Net Algorithm with Lazy Reduction on twsited KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_tw_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_tw_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_tw_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm on twsited KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_cz_tw_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_tw_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_cz_tw_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on twsited KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_cz_tw_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_tw_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_cz_tw_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv) on twsited KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_tw_ri_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_tw_ri_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_tw_ri_kss(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction on twsited KSS curves") {
		ep_rand(P);
	   	ep3_rand(Q);
		opt_en_tw_ri_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_tw_ri_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep3_set_infty(Q);
		opt_en_tw_ri_kss_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	printf("\nTesting if optimal ate pairing is bilinear:\n\n");

	TEST_BEGIN("The Elliptic Net Algorithm on twsited KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_tw_kss(F1,R,P,A,B,x,1);
		opt_en_tw_kss(F2,Q,P,A,B,x,1);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_tw_kss(F2,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The Elliptic Net Algorithm with Lazy Reduction on twsited KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_tw_kss_lz(F1,R,P,A,B,x,2);
		opt_en_tw_kss_lz(F2,Q,P,A,B,x,2);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_tw_kss_lz(F2,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm on twsited KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_cz_tw_kss(F1,R,P,A,B,x,1);
		opt_en_cz_tw_kss(F2,Q,P,A,B,x,1);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_tw_kss(F2,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on twsited KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_cz_tw_kss_lz(F1,R,P,A,B,x,2);
		opt_en_cz_tw_kss_lz(F2,Q,P,A,B,x,2);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_tw_kss_lz(F2,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv) on twsited KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_tw_ri_kss(F1,R,P,A,B,x,1);
		opt_en_tw_ri_kss(F2,Q,P,A,B,x,1);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_tw_ri_kss(F2,Q,P,A,B,x,1);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction on twsited KSS curves") {
		ep_rand(P);
		ep3_rand(Q);
		bn_rand_mod(s, r);
		ep3_mul_basic(R,Q,s);
		opt_en_tw_ri_kss_lz(F1,R,P,A,B,x,2);
		opt_en_tw_ri_kss_lz(F2,Q,P,A,B,x,2);
		fp18_exp(F2, F2, s);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_tw_ri_kss_lz(F2,Q,P,A,B,x,2);
		TEST_ASSERT(fp18_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
	end:

	fp_free(B1);
	fp3_free(A);
	fp3_free(B);
	bn_free(x);
	bn_free(r);
	bn_free(s);
	fp18_free(F1);
	fp18_free(F2);
	ep_free(P);
	ep3_free(Q);
	ep3_free(R);
}

static void pairing18(){
    fp_t B1;
	fp3_t A,B;
	fp18_t F1,F2;
	bn_t r,s,x;
	ep_t P;
	ep3_t Q,R;

  	fp_null(B1);
	fp3_null(A);
	fp3_null(B);
	bn_null(x);
	bn_null(r);
	bn_null(s);
	fp18_null(F1);
	fp18_null(F2);
	ep_null(P);
	ep3_null(Q);
	ep3_null(R);

	fp_new(B1);
	fp3_new(A);
	fp3_new(B);
	bn_new(x);
	bn_new(r);
	bn_new(s);
	fp18_new(F1);
	fp18_new(F2);
	ep_new(P);
	ep3_new(Q);
	ep3_new(R);

	ep_curve_get_ord(r);
	fp_prime_get_par(x);
	ep3_curve_get_a(A);//A=0
	ep3_curve_get_b(B);//B=-u^2
	fp_set_dig(B1,2);

/*==========================================================================*/
/*                       Easy efficency comparison                          */
/*==========================================================================*/

	BENCH_BEGIN("Miller's Algorithm") {
		ep_rand(P);
		ep3_rand(Q);
		BENCH_ADD(pp_map_oatep_k18(F1, P, Q));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm on KSS curves") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_kss(F1,Q,P,A[0],B1,x,1));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm with lazy reduction on KSS curves") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_kss_lz(F1,Q,P,A[0],B1,x,2));
	}
	BENCH_END;


	BENCH_BEGIN("The improved Elliptic Net Algorithm on KSS curves") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_cz_kss(F1,Q,P,A[0],B1,x,3));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on KSS curves") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_cz_kss_lz(F1,Q,P,A[0],B1,x,4));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  on KSS curves") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_cz_ri_kss(F1,Q,P,A[0],B1,x,1));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction KSS curves") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_cz_ri_kss_lz(F1,Q,P,A[0],B1,x,2));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm on twisted KSS curves") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_tw_kss(F1,Q,P,A,B,x,1));
	}
	BENCH_END;


	BENCH_BEGIN("The Elliptic Net Algorithm on twisted KSS curves with Lazy Reduction costs:") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_tw_kss_lz(F1,Q,P,A,B,x,2));
	}
	BENCH_END;


	BENCH_BEGIN("The improved Elliptic Net Algorithm on twisted KSS curves") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_cz_tw_kss(F1,Q,P,A,B,x,1));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm on twisted KSS curves with Lazy Reduction costs:") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_cz_tw_kss_lz(F1,Q,P,A,B,x,2));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm on twisted KSS curves (Eliminate Inv)  costs:") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_tw_ri_kss(F1,Q,P,A,B,x,1));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm on twisted KSS curves (Eliminate Inv)  with Lazy Reduction costs:") {
		ep_rand(P);
 		ep3_rand(Q);
		BENCH_ADD(opt_en_tw_ri_kss_lz(F1,Q,P,A,B,x,2));
	}
	BENCH_END;

	fp_free(B1);
	fp3_free(A);
	fp3_free(B);
	bn_free(x);
	bn_free(r);
	bn_free(s);
	fp18_free(F1);
	fp18_free(F2);
	ep_free(P);
	ep3_free(Q);
	ep3_free(R);
}

int main(){
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	conf_print();

	util_banner("Benchmarks for the PP module:", 0);

	if (ep_param_set_any_endom()!= RLC_OK) {
		THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

	ep_param_print();
	util_banner("Arithmetic:", 1);

	ep3_curve_init();
	ep3_curve_set_twist(1);//only can type (EP_DTYPE)1 or (EP_MTYPE)2,cof=1

	if (ep_param_embed() == 18) {
		if (pairing() == RLC_OK) {
			util_banner("All tests have passed.\n", 0);
			pairing18();
		}
		else
		{		
			core_clean();
			return 1;
		}
	}

	core_clean();
	return 0;
}
