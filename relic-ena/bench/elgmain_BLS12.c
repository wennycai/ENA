#include"ena.h"
#include"time.h"
//rdc rule:for add and sub, the two element should be >p or <p both;for mul and sqr,the two element should be <p both
/*BLS12-P381 curve parameter:
//y^2 = x^3 + 4
A=0;
B=4;
bata=-1;
Fq^2=Fq[u]/(u^2-beta);
ksai=u+1;
Fq^6=Fq^2[v]/(v^3-ksai);
Fq^12=Fp^6[w]/(w^2-v);
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
 * Benchmarks for the optimal ate pairing defined over the BLS12-P381 curve (elliptic net algorithm).
 *
 * @ingroup bench
 */

static int pairing(){
	int code = RLC_ERR;
    	fp_t B1;
	fp2_t A,B;
	fp12_t F1,F2;
	bn_t r,s,x;
	ep_t P;
	ep2_t Q,R;

  	fp_null(B1);
	fp2_null(A);
	fp2_null(B);
	bn_null(x);
	bn_null(r);
	bn_null(s);
	fp12_null(F1);
	fp12_null(F2);
	ep_null(P);
	ep2_null(Q);
	ep2_null(R);

	TRY{
	fp_new(B1);
	fp2_new(A);
	fp2_new(B);
	bn_new(x);
	bn_new(r);
	bn_new(s);
	fp12_new(F1);
	fp12_new(F2);
	ep_new(P);
	ep2_new(Q);
	ep2_new(R);

	ep_curve_get_ord(r);
	fp_prime_get_par(x);
	ep2_curve_get_a(A);//A=0
	ep2_curve_get_b(B);//B=4+4i
	fp_set_dig(B1,4);

/*==========================================================================*/
/*                                  Testing                                 */
/*==========================================================================*/
	printf("The correctness of the elliptic net algorithm on the untwist BLS curve:\n");

	printf("\nTesting if optimal ate pairing non-degeneracy is correct:\n\n");

	TEST_BEGIN("The Elliptic Net Algorithm on BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_bls(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_bls(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_bls(F1,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The Elliptic Net Algorithm with Lazy Reduction on BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_bls_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_bls_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_bls_lz(F1,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm on BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_cz_bls(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_bls(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_cz_bls(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_cz_bls_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_bls_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_cz_bls_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv) on BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_cz_ri_bls(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_ri_bls(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_cz_ri_bls(F1,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction on BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_cz_ri_bls_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_ri_bls_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_cz_ri_bls_lz(F1,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	printf("\nTesting if optimal ate pairing is bilinear:\n\n");

	TEST_BEGIN("The Elliptic Net Algorithm on BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_bls(F1,R,P,A[0],B1,x,1);
		opt_en_bls(F2,Q,P,A[0],B1,x,1);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_bls(F2,Q,P,A[0],B1,x,1);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The Elliptic Net Algorithm with Lazy Reduction on BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_bls_lz(F1,R,P,A[0],B1,x,2);
		opt_en_bls_lz(F2,Q,P,A[0],B1,x,2);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_bls_lz(F2,Q,P,A[0],B1,x,2);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm on BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_cz_bls(F1,R,P,A[0],B1,x,3);
		opt_en_cz_bls(F2,Q,P,A[0],B1,x,3);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_bls(F2,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_cz_bls_lz(F1,R,P,A[0],B1,x,4);
		opt_en_cz_bls_lz(F2,Q,P,A[0],B1,x,4);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_bls_lz(F2,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv) on BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_cz_ri_bls(F1,R,P,A[0],B1,x,3);
		opt_en_cz_ri_bls(F2,Q,P,A[0],B1,x,3);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_ri_bls(F2,Q,P,A[0],B1,x,3);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction on BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_cz_ri_bls_lz(F1,R,P,A[0],B1,x,4);
		opt_en_cz_ri_bls_lz(F2,Q,P,A[0],B1,x,4);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_ri_bls_lz(F2,Q,P,A[0],B1,x,4);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	printf("The correctness of the elliptic net algorithm on the twist BLS curve:\n");

	printf("\nTesting if optimal ate pairing non-degeneracy is correct:\n\n");

	TEST_BEGIN("The Elliptic Net Algorithm on twsited BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_tw_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_tw_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_tw_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The Elliptic Net Algorithm with Lazy Reduction on twsited BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_tw_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_tw_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_tw_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm on twsited BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_cz_tw_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_tw_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_cz_tw_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on twsited BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_cz_tw_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_cz_tw_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_cz_tw_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv) on twsited BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_tw_ri_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_tw_ri_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_tw_ri_bls(F1,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction on twsited BLS curves is non-degenerate") {
		ep_rand(P);
	   	ep2_rand(Q);
		opt_en_tw_ri_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) != RLC_EQ, end);
		ep_set_infty(P);
		opt_en_tw_ri_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
		ep_rand(P);
		ep2_set_infty(Q);
		opt_en_tw_ri_bls_lz(F1,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp_dig(F1, 1) == RLC_EQ, end);
	} TEST_END;

	printf("\nTesting if optimal ate pairing is bilinear:\n\n");

	TEST_BEGIN("The Elliptic Net Algorithm on twsited BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_tw_bls(F1,R,P,A,B,x,1);
		opt_en_tw_bls(F2,Q,P,A,B,x,1);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_tw_bls(F2,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The Elliptic Net Algorithm with Lazy Reduction on twsited BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_tw_bls_lz(F1,R,P,A,B,x,2);
		opt_en_tw_bls_lz(F2,Q,P,A,B,x,2);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_tw_bls_lz(F2,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm on twsited BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_cz_tw_bls(F1,R,P,A,B,x,1);
		opt_en_cz_tw_bls(F2,Q,P,A,B,x,1);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_tw_bls(F2,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on twsited BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_cz_tw_bls_lz(F1,R,P,A,B,x,2);
		opt_en_cz_tw_bls_lz(F2,Q,P,A,B,x,2);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_cz_tw_bls_lz(F2,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv) on twsited BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_tw_ri_bls(F1,R,P,A,B,x,1);
		opt_en_tw_ri_bls(F2,Q,P,A,B,x,1);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_tw_ri_bls(F2,Q,P,A,B,x,1);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;

	TEST_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction on twsited BLS curves is bilinear") {
		ep_rand(P);
		ep2_rand(Q);
		bn_rand_mod(s, r);
		ep2_mul(R,Q,s);
		opt_en_tw_ri_bls_lz(F1,R,P,A,B,x,2);
		opt_en_tw_ri_bls_lz(F2,Q,P,A,B,x,2);
		fp12_exp(F2, F2, s);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
		ep_mul(P, P, s);
		opt_en_tw_ri_bls_lz(F2,Q,P,A,B,x,2);
		TEST_ASSERT(fp12_cmp(F1, F2) == RLC_EQ, end);
	} TEST_END;
	}
		CATCH_ANY{
			util_print("FATAL ERROR!\n");
			ERROR(end);
	}
	code = RLC_OK;
	end:

	fp_free(B1);
	fp2_free(A);
	fp2_free(B);
	bn_free(x);
	bn_free(r);
	bn_free(s);
	fp12_free(F1);
	fp12_free(F2);
	ep_free(P);
	ep2_free(Q);
	ep2_free(R);
}

static void pairing12(){
   	fp_t B1;
	fp2_t A,B;
	fp12_t F1,F2;
	bn_t r,s,x;
	ep_t P;
	ep2_t Q,R;

  	fp_null(B1);
	fp2_null(A);
	fp2_null(B);
	bn_null(x);
	bn_null(r);
	bn_null(s);
	fp12_null(F1);
	fp12_null(F2);
	ep_null(P);
	ep2_null(Q);
	ep2_null(R);

	fp_new(B1);
	fp2_new(A);
	fp2_new(B);
	bn_new(x);
	bn_new(r);
	bn_new(s);
	fp12_new(F1);
	fp12_new(F2);
	ep_new(P);
	ep2_new(Q);
	ep2_new(R);

	ep_curve_get_ord(r);
	fp_prime_get_par(x);
	ep2_curve_get_a(A);//A=0
	ep2_curve_get_b(B);//B=4+4i
	fp_set_dig(B1,4);

/*==========================================================================*/
/*                       Easy efficency comparison                          */
/*==========================================================================*/

	BENCH_BEGIN("Miller's Algorithm") {
		ep_rand(P);
		ep2_rand(Q);
		BENCH_ADD(pp_map_oatep_k12(F1, P, Q));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm on BLS curves") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_bls(F1,Q,P,A[0],B1,x,1));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm with lazy reduction on BLS curves") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_bls_lz(F1,Q,P,A[0],B1,x,2));
	}
	BENCH_END;


	BENCH_BEGIN("The improved Elliptic Net Algorithm on BLS curves") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_cz_bls(F1,Q,P,A[0],B1,x,3));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm with Lazy Reduction on BLS curves") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_cz_bls_lz(F1,Q,P,A[0],B1,x,4));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  on BLS curves") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_cz_ri_bls(F1,Q,P,A[0],B1,x,3));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm (Eliminate Inv)  with Lazy Reduction BLS curves") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_cz_ri_bls_lz(F1,Q,P,A[0],B1,x,4));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm on twisted BLS curves") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_tw_bls(F1,Q,P,A,B,x,1));
	}
	BENCH_END;


	BENCH_BEGIN("The Elliptic Net Algorithm on twisted BLS curves with Lazy Reduction costs:") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_tw_bls_lz(F1,Q,P,A,B,x,2));
	}
	BENCH_END;


	BENCH_BEGIN("The improved Elliptic Net Algorithm on twisted BLS curves") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_cz_tw_bls(F1,Q,P,A,B,x,1));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm on twisted BLS curves with Lazy Reduction costs:") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_cz_tw_bls_lz(F1,Q,P,A,B,x,2));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm on twisted BLS curves (Eliminate Inv)  costs:") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_tw_ri_bls(F1,Q,P,A,B,x,1));
	}
	BENCH_END;

	BENCH_BEGIN("The Elliptic Net Algorithm on twisted BLS curves (Eliminate Inv)  with Lazy Reduction costs:") {
		ep_rand(P);
 		ep2_rand(Q);
		BENCH_ADD(opt_en_tw_ri_bls_lz(F1,Q,P,A,B,x,2));
	}
	BENCH_END;

	fp_free(B1);
	fp2_free(A);
	fp2_free(B);
	bn_free(x);
	bn_free(r);
	bn_free(s);
	fp12_free(F1);
	fp12_free(F2);
	ep_free(P);
	ep2_free(Q);
	ep2_free(R);
}


int main(){
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	conf_print();

	util_banner("Benchmarks for the PP module:", 0);

	if (ep_param_set_any_pairf() != RLC_OK) {
		THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

	ep_param_print();
	util_banner("Arithmetic:", 1);

	if (ep_param_embed() == 12) {
		if (pairing() == RLC_OK) {
			util_banner("All tests have passed.\n", 0);
			pairing12();
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
