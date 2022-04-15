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
 * Benchmarks for scalar multiplication defined over the NIST-P384 curve (elliptic net algorithm).
 *
 * @ingroup bench
 */

/* We chose a base point P (generator) */

#define NIST_P384_A	"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFC"
#define NIST_P384_B	"B3312FA7E23EE7E4988E056BE3F82D19181D9C6EFE8141120314088F5013875AC656398D8A2ED19D2A85C8EDD3EC2AEF"
#define xp "61784D9786939510C2F522AD4925E08C4C23BAD511BAFC05BC9C295B2C61BC009E08355AF9785AE265ABE32253CCE4F8"
#define yp "7C64BF118E32A13C14DE4A26542FF87C13F4CE47F0157307027E1974D6E25AA64AA025671B1317970DCA68C318EC5E59"
#define cuberoot "2259937780625989193807398058089059511128447613510482027874020612867561784561402210729349239648684610148283527454122"

static int scalmul(){
	int code = RLC_ERR;
	fp_t A,B,C;	
	bn_t r,s;
	ep_t P,mP,sP;

	fp_null(A);
	fp_null(B);
	fp_null(C);
	bn_null(r);
	bn_null(s);
	ep_null(P);
	ep_null(mP);
	ep_null(sP);

	TRY{
		fp_new(A);
		fp_new(B);
		fp_new(C);
		bn_new(r);
		bn_new(s);
		ep_new(P);
		ep_new(mP);
		ep_new(sP);		

		fp_read_str(A,NIST_P384_A,strlen(NIST_P384_A),16);
		fp_read_str(B,NIST_P384_B,strlen(NIST_P384_B),16);

		ep_curve_get_ord(r);
		ep_rand(P);

		fp_read_str(P->x,xp,strlen(xp),16);
		fp_read_str(P->y,yp,strlen(yp),16);	
		fp_read_str(C,cuberoot,strlen(cuberoot),10);

/*==========================================================================*/
/*                                  Testing                                 */
/*==========================================================================*/
	
		printf("\nTesting if scalar multiplication is correct \n");

		TEST_BEGIN("The improved Elliptic Net Algorithm for scalar multiplication is correct") {
			bn_rand_mod(s, r);
			ep_mul(sP,P,s);
			scalmul_net_cz(mP,P,A,B,s,C,1);
			TEST_ASSERT(ep_cmp(sP, mP) == RLC_EQ, end);
			scalmul_net_cz(mP,P,A,B,s,C,2);
			TEST_ASSERT(ep_cmp(sP, mP) == RLC_EQ, end);
		} TEST_END;

		TEST_BEGIN("The updated Elliptic Net Algorithm for scalar multiplication is correct") {
			bn_rand_mod(s, r);
			ep_mul(sP,P,s);
			scalmul_net_cz(mP,P,A,B,s,C,3);
			TEST_ASSERT(ep_cmp(sP, mP) == RLC_EQ, end);
			scalmul_net_cz(mP,P,A,B,s,C,4);
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
	fp_free(C);
	bn_free(r);
	bn_free(s);
	ep_free(P);
	ep_free(mP);
	ep_free(sP);		
}

static void multiplication(){
	fp_t A,B,C;	
	bn_t r,s;
	ep_t P,mP,sP;

	fp_null(A);
	fp_null(B);
	fp_null(C);
	bn_null(r);
	bn_null(s);
	ep_null(P);
	ep_null(mP);
	ep_null(sP);

	fp_new(A);
	fp_new(B);
	fp_new(C);
	bn_new(r);
	bn_new(s);
	ep_new(P);
	ep_new(mP);
	ep_new(sP);		

	fp_read_str(A,NIST_P384_A,strlen(NIST_P384_A),16);
	fp_read_str(B,NIST_P384_B,strlen(NIST_P384_B),16);

	ep_curve_get_ord(r);
		
	ep_rand(P);

	fp_read_str(P->x,xp,strlen(xp),16);
	fp_read_str(P->y,yp,strlen(yp),16);	
	fp_read_str(C,cuberoot,strlen(cuberoot),10);

/*==========================================================================*/
/*                       Easy efficency comparison                          */
/*==========================================================================*/

	BENCH_BEGIN("The improved Elliptic Net Algorithm for scalar multiplication cost") {
		bn_rand_mod(s, r);
		BENCH_ADD(scalmul_net_cz(sP,P,A,B,s,C,1));
	}
	BENCH_END;

	BENCH_BEGIN("The improved Elliptic Net Algorithm for scalar multiplication with lazy reduction cost") {
		bn_rand_mod(s, r);
		BENCH_ADD(scalmul_net_cz(sP,P,A,B,s,C,2));
	}
	BENCH_END;

	BENCH_BEGIN("The updated Elliptic Net Algorithm for scalar multiplication cost") {
		bn_rand_mod(s, r);
		BENCH_ADD(scalmul_net_cz(sP,P,A,B,s,C,3));
	}
	BENCH_END;

	BENCH_BEGIN("The updated Elliptic Net Algorithm for scalar multiplication with lazy reduction cost") {
		bn_rand_mod(s, r);
		BENCH_ADD(scalmul_net_cz(sP,P,A,B,s,C,4));
	}
	BENCH_END;

	fp_free(A);
	fp_free(B);
	fp_free(C);
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

	ep_param_set(NIST_P384);

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
