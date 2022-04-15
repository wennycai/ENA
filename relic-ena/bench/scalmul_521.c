#include"ena.h"
#include"time.h"

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
 * Benchmarks for scalar multiplication defined over the NIST-P521 curve (elliptic net algorithm).
 *
 * @ingroup bench
 */

#define NIST_P521_A	"1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC"
#define NIST_P521_B	"51953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00"

static int scalmul(){
	int code = RLC_ERR;
	fp_t A,B;	
	bn_t r,s;
	ep_t P,mP,sP;

	fp_null(A);
	fp_null(B);
	bn_null(r);
	bn_null(s);
	ep_null(P);
	ep_null(mP);
	ep_null(sP);

	TRY{
		fp_new(A);
		fp_new(B);
		bn_new(r);
		bn_new(s);
		ep_new(P);
		ep_new(mP);
		ep_new(sP);		

		fp_read_str(A,NIST_P521_A,strlen(NIST_P521_A),16);
		fp_read_str(B,NIST_P521_B,strlen(NIST_P521_B),16);

		ep_curve_get_ord(r);
		ep_rand(P);

/*==========================================================================*/
/*                                  Testing                                 */
/*==========================================================================*/
	
		printf("\nTesting if scalar multiplication is correct \n");

		TEST_BEGIN("The improved Elliptic Net Algorithm for scalar multiplication is correct") {
			bn_rand_mod(s, r);
			ep_mul(sP,P,s);
			scalmul_net(mP,P,A,B,s,1);
			TEST_ASSERT(ep_cmp(sP, mP) == RLC_EQ, end);
			scalmul_net(mP,P,A,B,s,2);
			TEST_ASSERT(ep_cmp(sP, mP) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("The updated Elliptic Net Algorithm for scalar multiplication is correct") {
			bn_rand_mod(s, r);
			ep_mul(sP,P,s);
			scalmul_net(mP,P,A,B,s,3);
			TEST_ASSERT(ep_cmp(sP, mP) == RLC_EQ, end);
			scalmul_net(mP,P,A,B,s,4);
			TEST_ASSERT(ep_cmp(sP, mP) == RLC_EQ, end);
		} TEST_END;
	}
	CATCH_ANY{
		util_print("FATAL ERROR!\n");
		ERROR(end);
	}
	code = RLC_OK;
	end:

	fp_free(A);
	fp_free(B);
	bn_free(r);
	bn_free(s);
	ep_free(P);
	ep_free(mP);
	ep_free(sP);		
}

static void multiplication(){
	fp_t A,B;	
	bn_t r,s;
	ep_t P,mP,sP;

	fp_null(A);
	fp_null(B);
	bn_null(r);
	bn_null(s);
	ep_null(P);
	ep_null(mP);
	ep_null(sP);

	fp_new(A);
	fp_new(B);
	bn_new(r);
	bn_new(s);
	ep_new(P);
	ep_new(mP);
	ep_new(sP);		

	fp_read_str(A,NIST_P521_A,strlen(NIST_P521_A),16);
	fp_read_str(B,NIST_P521_B,strlen(NIST_P521_B),16);

	ep_curve_get_ord(r);
		
	ep_rand(P);

/*==========================================================================*/
/*                       Easy efficency comparison                          */
/*==========================================================================*/

	BENCH_BEGIN("The improved Elliptic Net Algorithm for scalar multiplication cost") {
		bn_rand_mod(s, r);
		BENCH_ADD(scalmul_net(sP,P,A,B,s,1));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm for scalar multiplication with lazy reduction cost") {
		bn_rand_mod(s, r);
		BENCH_ADD(scalmul_net(sP,P,A,B,s,2));
	}
	BENCH_END;

	BENCH_BEGIN("The updated Elliptic Net Algorithm for scalar multiplication cost") {
		bn_rand_mod(s, r);
		BENCH_ADD(scalmul_net(sP,P,A,B,s,3));
	}
	BENCH_END;

	BENCH_BEGIN("The updated Elliptic Net Algorithm for scalar multiplication with lazy reduction cost") {
		bn_rand_mod(s, r);
		BENCH_ADD(scalmul_net(sP,P,A,B,s,4));
	}
	BENCH_END;

	fp_free(A);
	fp_free(B);
	bn_free(r);
	bn_free(s);
	ep_free(P);
	ep_free(mP);
	ep_free(sP);	
}


int main(){
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	conf_print();

	util_banner("Benchmarks for the PP module:", 0);

	ep_param_set(NIST_P521);

	ep_param_print();

	util_banner("Arithmetic:", 1);

	if (scalmul() == RLC_OK) {
		util_banner("All tests have passed.\n", 0);
		multiplication();
	}
	else
	{		
		core_clean();
		return 1;
	}

	core_clean();
	return 0;
}
