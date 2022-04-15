/*==========================================================================*/
/*===========                                              =================*/
/*            Scalar MUltiplication Based on Elliptic Nets                  */
/*===========                                              =================*/
/*==========================================================================*/

#include"ena.h"


/*improved double and add alg by Rao et al.*/
static void DoubleAndAdd01(fp_t V1[],fp_t S1[],fp_t P1[],fp_t k[],fp_t alpha,bn_t m){
	int nb,i;	
	nb = bn_bits(m);
	nb = nb - 1;
	
	for(i=nb-1;i>=0;i--){
		int b=bn_get_bit(m,i);	
		for (int j = 0; j <= 5; j++)//5S,5M
		{
			fp_sqr(S1[j],V1[j + 1]);
			fp_mul(P1[j],V1[j],V1[j + 2]);
		}
		
		if (b == 0)
		{
			fp_add(k[0],P1[0],P1[1]);
			fp_sub(k[5],S1[0],S1[1]);
			fp_mul(k[0],k[0],k[5]);

			fp_sub(k[1],P1[0],P1[1]);
			fp_add(k[5],S1[0],S1[1]);
			fp_mul(k[1],k[1],k[5]);

			fp_add(k[2],P1[0],P1[2]);
			fp_sub(k[5],S1[0],S1[2]);
			fp_mul(k[2],k[2],k[5]);

			fp_sub(k[3],P1[0],P1[2]);
			fp_add(k[5],S1[0],S1[2]);
			fp_mul(k[3],k[3],k[5]);

			fp_sub(k[4],P1[1],P1[2]);
			fp_add(k[5],S1[1],S1[2]);
			fp_mul(k[4],k[4],k[5]);
			fp_dbl(k[4],k[4]);

			fp_sub(V1[0],k[0],k[1]);
			fp_hlv_integ(V1[0],V1[0]);
			fp_sub(V1[1],k[2],k[3]);
			fp_hlv_integ(V1[1],V1[1]);
			fp_mul(V1[1],V1[1],alpha);

			fp_add(V1[2],k[2],k[3]);
			fp_add(k[5],k[0],k[1]);
			fp_add(k[5],k[5],k[4]);
			fp_sub(V1[2],V1[2],k[5]);
			fp_hlv_integ(V1[2],V1[2]);

			fp_mul(V1[3],S1[1],P1[3]);
			fp_mul(k[5],S1[3],P1[1]);
			fp_sub(V1[3],V1[3],k[5]);
			fp_mul(V1[3],V1[3],alpha);


			fp_add(k[0],P1[2],P1[3]);
			fp_sub(k[5],S1[2],S1[3]);
			fp_mul(k[0],k[0],k[5]);

			fp_sub(k[1],P1[2],P1[3]);
			fp_add(k[5],S1[2],S1[3]);
			fp_mul(k[1],k[1],k[5]);

			fp_add(k[2],P1[2],P1[4]);
			fp_sub(k[5],S1[2],S1[4]);
			fp_mul(k[2],k[2],k[5]);

			fp_sub(k[3],P1[2],P1[4]);
			fp_add(k[5],S1[2],S1[4]);
			fp_mul(k[3],k[3],k[5]);

			fp_sub(k[4],P1[3],P1[4]);
			fp_add(k[5],S1[3],S1[4]);
			fp_mul(k[4],k[4],k[5]);
			fp_dbl(k[4],k[4]);

			fp_sub(V1[4],k[0],k[1]);
			fp_hlv_integ(V1[4],V1[4]);
			fp_sub(V1[5],k[2],k[3]);
			fp_hlv_integ(V1[5],V1[5]);
			fp_mul(V1[5],V1[5],alpha);
			fp_add(V1[6],k[2],k[3]);
			fp_add(k[5],k[0],k[1]);
			fp_add(k[5],k[5],k[4]);
			fp_sub(V1[6],V1[6],k[5]);
			fp_hlv_integ(V1[6],V1[6]);
			fp_mul(V1[7],S1[3],P1[5]);
			fp_mul(k[5],S1[5],P1[3]);
			fp_sub(V1[7],V1[7],k[5]);
			fp_mul(V1[7],V1[7],alpha);

		}
		else
		{
			fp_add(k[0],P1[1],P1[2]);
			fp_sub(k[5],S1[1],S1[2]);
			fp_mul(k[0],k[0],k[5]);

			fp_sub(k[1],P1[1],P1[2]);
			fp_add(k[5],S1[1],S1[2]);
			fp_mul(k[1],k[1],k[5]);

			fp_add(k[2],P1[1],P1[3]);
			fp_sub(k[5],S1[1],S1[3]);
			fp_mul(k[2],k[2],k[5]);

			fp_sub(k[3],P1[1],P1[3]);
			fp_add(k[5],S1[1],S1[3]);
			fp_mul(k[3],k[3],k[5]);

			fp_sub(k[4],P1[2],P1[3]);
			fp_add(k[5],S1[2],S1[3]);
			fp_mul(k[4],k[4],k[5]);
			fp_dbl(k[4],k[4]);

			fp_mul(V1[0],S1[0],P1[2]);
			fp_mul(k[5],S1[2],P1[0]);
			fp_sub(V1[0],V1[0],k[5]);
			fp_mul(V1[0],V1[0],alpha);
			fp_sub(V1[1],k[0],k[1]);
			fp_hlv_integ(V1[1],V1[1]);
			fp_sub(V1[2],k[2],k[3]);
			fp_hlv_integ(V1[2],V1[2]);
			fp_mul(V1[2],V1[2],alpha);
			fp_add(V1[3],k[2],k[3]);
			fp_add(k[5],k[0],k[1]);
			fp_add(k[5],k[5],k[4]);
			fp_sub(V1[3],V1[3],k[5]);
			fp_hlv_integ(V1[3],V1[3]);
			fp_mul(V1[4],S1[2],P1[4]);
			fp_mul(k[5],S1[4],P1[2]);
			fp_sub(V1[4],V1[4],k[5]);
			fp_mul(V1[4],V1[4],alpha);

			fp_add(k[0],P1[3],P1[4]);
			fp_sub(k[5],S1[3],S1[4]);
			fp_mul(k[0],k[0],k[5]);

			fp_sub(k[1],P1[3],P1[4]);
			fp_add(k[5],S1[3],S1[4]);
			fp_mul(k[1],k[1],k[5]);

			fp_add(k[2],P1[3],P1[5]);
			fp_sub(k[5],S1[3],S1[5]);
			fp_mul(k[2],k[2],k[5]);

			fp_sub(k[3],P1[3],P1[5]);
			fp_add(k[5],S1[3],S1[5]);
			fp_mul(k[3],k[3],k[5]);

			fp_sub(k[4],P1[4],P1[5]);
			fp_add(k[5],S1[4],S1[5]);
			fp_mul(k[4],k[4],k[5]);
			fp_dbl(k[4],k[4]);

			fp_sub(V1[5],k[0],k[1]);
			fp_hlv_integ(V1[5],V1[5]);
			fp_sub(V1[6],k[2],k[3]);
			fp_hlv_integ(V1[6],V1[6]);
			fp_mul(V1[6],V1[6],alpha);
			fp_add(V1[7],k[2],k[3]);
			fp_add(k[5],k[0],k[1]);
			fp_add(k[5],k[5],k[4]);
			fp_sub(V1[7],V1[7],k[5]);
			fp_hlv_integ(V1[7],V1[7]);
		}
	}
}

/*improved double and add alg by Rao et al.*/
static void DoubleAndAdd01_lz(fp_t V1[],fp_t S1[],fp_t P1[],fp_t k[],fp_t alpha,bn_t m){
	int nb,i;	
	nb = bn_bits(m);
	nb = nb - 1;
	dv_t t[4],T0;

	dv_null(T0);
	for(i=0;i<4;i++)dv_null(t[i]);
	
	dv_new(T0);
	for(i=0;i<4;i++)dv_new(t[i]);
	
	for(i=nb-1;i>=0;i--){
		int b=bn_get_bit(m,i);	
		for (int j = 0; j <= 5; j++)//5S,5M
		{
			fp_sqr(S1[j],V1[j + 1]);
			fp_mul(P1[j],V1[j],V1[j + 2]);
		}
		
		if (b == 0)
		{
			fp_add(k[0],P1[0],P1[1]);
			fp_sub(k[5],S1[0],S1[1]);
			fp_muln_low(t[0],k[0],k[5]);

			fp_sub(k[1],P1[0],P1[1]);
			fp_add(k[5],S1[0],S1[1]);
			fp_muln_low(t[1],k[1],k[5]);

			fp_add(k[2],P1[0],P1[2]);
			fp_sub(k[5],S1[0],S1[2]);
			fp_muln_low(t[2],k[2],k[5]);

			fp_sub(k[3],P1[0],P1[2]);
			fp_add(k[5],S1[0],S1[2]);
			fp_muln_low(t[3],k[3],k[5]);

			fp_sub(k[4],P1[1],P1[2]);
			fp_add(k[5],S1[1],S1[2]);
			fp_mul(k[4],k[4],k[5]);

			fp_subc_low(T0,t[0],t[1]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[0],T0);
			fp_subc_low(T0,t[2],t[3]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[1],T0);
			fp_mul(V1[1],V1[1],alpha);


			fp_addc_low(T0,t[2],t[3]);
			fp_addc_low(t[2],t[0],t[1]);
			fp_subc_low(T0,T0,t[2]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[2],T0);
			fp_sub(V1[2],V1[2],k[4]);

			fp_muln_low(t[0],S1[1],P1[3]);
			fp_muln_low(t[1],S1[3],P1[1]);
			fp_subc_low(T0,t[0],t[1]);
			fp_rdc(V1[3],T0);
			fp_mul(V1[3],V1[3],alpha);


			fp_add(k[0],P1[2],P1[3]);
			fp_sub(k[5],S1[2],S1[3]);
			fp_muln_low(t[0],k[0],k[5]);

			fp_sub(k[1],P1[2],P1[3]);
			fp_add(k[5],S1[2],S1[3]);
			fp_muln_low(t[1],k[1],k[5]);

			fp_add(k[2],P1[2],P1[4]);
			fp_sub(k[5],S1[2],S1[4]);
			fp_muln_low(t[2],k[2],k[5]);

			fp_sub(k[3],P1[2],P1[4]);
			fp_add(k[5],S1[2],S1[4]);
			fp_muln_low(t[3],k[3],k[5]);

			fp_sub(k[4],P1[3],P1[4]);
			fp_add(k[5],S1[3],S1[4]);
			fp_mul(k[4],k[4],k[5]);
			//fp_dbl(k[4],k[4]);

			fp_subc_low(T0,t[0],t[1]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[4],T0);
			fp_subc_low(T0,t[2],t[3]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[5],T0);
			fp_mul(V1[5],V1[5],alpha);
			fp_addc_low(T0,t[2],t[3]);
			fp_addc_low(t[2],t[0],t[1]);
			fp_subc_low(T0,T0,t[2]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[6],T0);
			fp_sub(V1[6],V1[6],k[4]);
			fp_muln_low(t[0],S1[3],P1[5]);
			fp_muln_low(t[1],S1[5],P1[3]);
			fp_subc_low(T0,t[0],t[1]);
			fp_rdc(V1[7],T0);
			fp_mul(V1[7],V1[7],alpha);

		}
		else
		{
			fp_muln_low(t[0],S1[0],P1[2]);
			fp_muln_low(t[1],S1[2],P1[0]);
			fp_subc_low(T0,t[0],t[1]);
			fp_rdc(V1[0],T0);
			fp_mul(V1[0],V1[0],alpha);

			fp_add(k[0],P1[1],P1[2]);
			fp_sub(k[5],S1[1],S1[2]);
			fp_muln_low(t[0],k[0],k[5]);

			fp_sub(k[1],P1[1],P1[2]);
			fp_add(k[5],S1[1],S1[2]);
			fp_muln_low(t[1],k[1],k[5]);

			fp_add(k[2],P1[1],P1[3]);
			fp_sub(k[5],S1[1],S1[3]);
			fp_muln_low(t[2],k[2],k[5]);

			fp_sub(k[3],P1[1],P1[3]);
			fp_add(k[5],S1[1],S1[3]);
			fp_muln_low(t[3],k[3],k[5]);

			fp_sub(k[4],P1[2],P1[3]);
			fp_add(k[5],S1[2],S1[3]);
			fp_mul(k[4],k[4],k[5]);
			//fp_dbl(k[4],k[4]);

			fp_subc_low(T0,t[0],t[1]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[1],T0);
			fp_subc_low(T0,t[2],t[3]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[2],T0);
			fp_mul(V1[2],V1[2],alpha);
			fp_addc_low(T0,t[2],t[3]);
			fp_addc_low(t[2],t[0],t[1]);
			fp_subc_low(T0,T0,t[2]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[3],T0);
			fp_sub(V1[3],V1[3],k[4]);
			fp_muln_low(t[0],S1[2],P1[4]);
			fp_muln_low(t[1],S1[4],P1[2]);
			fp_subc_low(T0,t[0],t[1]);
			fp_rdc(V1[4],T0);
			fp_mul(V1[4],V1[4],alpha);

			fp_add(k[0],P1[3],P1[4]);
			fp_sub(k[5],S1[3],S1[4]);
			fp_muln_low(t[0],k[0],k[5]);

			fp_sub(k[1],P1[3],P1[4]);
			fp_add(k[5],S1[3],S1[4]);
			fp_muln_low(t[1],k[1],k[5]);

			fp_add(k[2],P1[3],P1[5]);
			fp_sub(k[5],S1[3],S1[5]);
			fp_muln_low(t[2],k[2],k[5]);

			fp_sub(k[3],P1[3],P1[5]);
			fp_add(k[5],S1[3],S1[5]);
			fp_muln_low(t[3],k[3],k[5]);

			fp_sub(k[4],P1[4],P1[5]);
			fp_add(k[5],S1[4],S1[5]);
			fp_mul(k[4],k[4],k[5]);
			//fp_dbl(k[4],k[4]);

			fp_subc_low(T0,t[0],t[1]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[5],T0);
			fp_subc_low(T0,t[2],t[3]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[6],T0);
			fp_mul(V1[6],V1[6],alpha);
			fp_addc_low(T0,t[2],t[3]);
			fp_addc_low(t[2],t[0],t[1]);
			fp_subc_low(T0,T0,t[2]);
			fp_hlvd_low(T0,T0);
			fp_rdc(V1[7],T0);
			fp_sub(V1[7],V1[7],k[4]);
		}
	}
	
	dv_free(T0);
	for(i=0;i<4;i++)dv_free(t[i]);
}

/*improved double and add alg*/
static void DoubleAndAdd02(fp_t V1[],fp_t S1[],fp_t P1[],fp_t k[],fp_t alpha,bn_t m){
	int nb,i;	
	nb = bn_bits(m);
	nb = nb - 1;
	
	for(i=nb-1;i>=0;i--){
		int b=bn_get_bit(m,i);	
		for (int j = 0; j <= 5; j++)//5S,5M
		{
			fp_sqr(S1[j],V1[j + 1]);
			fp_mul(P1[j],V1[j],V1[j + 2]);
		}
		
		if (b == 0)
		{
			fp_add(k[0],P1[0],P1[1]);
			fp_sub(k[5],S1[0],S1[1]);
			fp_mul(k[0],k[0],k[5]);

			fp_sub(k[1],P1[0],P1[1]);
			fp_add(k[5],S1[0],S1[1]);
			fp_mul(k[1],k[1],k[5]);

			fp_add(k[2],P1[0],P1[2]);
			fp_sub(k[5],S1[0],S1[2]);
			fp_mul(k[2],k[2],k[5]);

			fp_sub(k[3],P1[0],P1[2]);
			fp_add(k[5],S1[0],S1[2]);
			fp_mul(k[3],k[3],k[5]);

			fp_sub(k[4],P1[1],P1[2]);
			fp_add(k[5],S1[1],S1[2]);
			fp_mul(k[4],k[4],k[5]);
			fp_dbl(k[4],k[4]);

			fp_sub(V1[0],k[0],k[1]);
			fp_sub(V1[1],k[2],k[3]);
			fp_mul(V1[1],V1[1],alpha);
			fp_add(V1[2],k[2],k[3]);
			fp_add(k[5],k[0],k[1]);
			fp_add(k[5],k[5],k[4]);
			fp_sub(V1[2],V1[2],k[5]);
			fp_mul(V1[3],S1[1],P1[3]);
			fp_mul(k[5],S1[3],P1[1]);
			fp_sub(V1[3],V1[3],k[5]);
			fp_mul(V1[3],V1[3],alpha);
			fp_dbl(V1[3],V1[3]);


			fp_add(k[0],P1[2],P1[3]);
			fp_sub(k[5],S1[2],S1[3]);
			fp_mul(k[0],k[0],k[5]);

			fp_sub(k[1],P1[2],P1[3]);
			fp_add(k[5],S1[2],S1[3]);
			fp_mul(k[1],k[1],k[5]);

			fp_add(k[2],P1[2],P1[4]);
			fp_sub(k[5],S1[2],S1[4]);
			fp_mul(k[2],k[2],k[5]);

			fp_sub(k[3],P1[2],P1[4]);
			fp_add(k[5],S1[2],S1[4]);
			fp_mul(k[3],k[3],k[5]);

			fp_sub(k[4],P1[3],P1[4]);
			fp_add(k[5],S1[3],S1[4]);
			fp_mul(k[4],k[4],k[5]);
			fp_dbl(k[4],k[4]);

			fp_sub(V1[4],k[0],k[1]);
			fp_sub(V1[5],k[2],k[3]);
			fp_mul(V1[5],V1[5],alpha);
			fp_add(V1[6],k[2],k[3]);
			fp_add(k[5],k[0],k[1]);
			fp_add(k[5],k[5],k[4]);
			fp_sub(V1[6],V1[6],k[5]);
			fp_mul(V1[7],S1[3],P1[5]);
			fp_mul(k[5],S1[5],P1[3]);
			fp_sub(V1[7],V1[7],k[5]);
			fp_mul(V1[7],V1[7],alpha);
			fp_dbl(V1[7],V1[7]);

		}
		else
		{
			fp_add(k[0],P1[1],P1[2]);
			fp_sub(k[5],S1[1],S1[2]);
			fp_mul(k[0],k[0],k[5]);

			fp_sub(k[1],P1[1],P1[2]);
			fp_add(k[5],S1[1],S1[2]);
			fp_mul(k[1],k[1],k[5]);

			fp_add(k[2],P1[1],P1[3]);
			fp_sub(k[5],S1[1],S1[3]);
			fp_mul(k[2],k[2],k[5]);

			fp_sub(k[3],P1[1],P1[3]);
			fp_add(k[5],S1[1],S1[3]);
			fp_mul(k[3],k[3],k[5]);

			fp_sub(k[4],P1[2],P1[3]);
			fp_add(k[5],S1[2],S1[3]);
			fp_mul(k[4],k[4],k[5]);
			fp_dbl(k[4],k[4]);

			fp_mul(V1[0],S1[0],P1[2]);
			fp_mul(k[5],S1[2],P1[0]);
			fp_sub(V1[0],V1[0],k[5]);
			fp_mul(V1[0],V1[0],alpha);
			fp_dbl(V1[0],V1[0]);
			fp_sub(V1[1],k[0],k[1]);
			fp_sub(V1[2],k[2],k[3]);
			fp_mul(V1[2],V1[2],alpha);
			fp_add(V1[3],k[2],k[3]);
			fp_add(k[5],k[0],k[1]);
			fp_add(k[5],k[5],k[4]);
			fp_sub(V1[3],V1[3],k[5]);
			fp_mul(V1[4],S1[2],P1[4]);
			fp_mul(k[5],S1[4],P1[2]);
			fp_sub(V1[4],V1[4],k[5]);
			fp_mul(V1[4],V1[4],alpha);
			fp_dbl(V1[4],V1[4]);


			fp_add(k[0],P1[3],P1[4]);
			fp_sub(k[5],S1[3],S1[4]);
			fp_mul(k[0],k[0],k[5]);

			fp_sub(k[1],P1[3],P1[4]);
			fp_add(k[5],S1[3],S1[4]);
			fp_mul(k[1],k[1],k[5]);

			fp_add(k[2],P1[3],P1[5]);
			fp_sub(k[5],S1[3],S1[5]);
			fp_mul(k[2],k[2],k[5]);

			fp_sub(k[3],P1[3],P1[5]);
			fp_add(k[5],S1[3],S1[5]);
			fp_mul(k[3],k[3],k[5]);

			fp_sub(k[4],P1[4],P1[5]);
			fp_add(k[5],S1[4],S1[5]);
			fp_mul(k[4],k[4],k[5]);
			fp_dbl(k[4],k[4]);

			fp_sub(V1[5],k[0],k[1]);
			fp_sub(V1[6],k[2],k[3]);
			fp_mul(V1[6],V1[6],alpha);
			fp_add(V1[7],k[2],k[3]);
			fp_add(k[5],k[0],k[1]);
			fp_add(k[5],k[5],k[4]);
			fp_sub(V1[7],V1[7],k[5]);
		}
	}
}

/*improved double and add alg with lazy reduction*/
static void DoubleAndAdd02_lz(fp_t V1[],fp_t S1[],fp_t P1[],fp_t k[],fp_t alpha,bn_t m){
	int nb,i;	
	nb = bn_bits(m);
	nb = nb - 1;
	dv_t t[5],T0;

	dv_null(T0);
	for(i=0;i<5;i++)dv_null(t[i]);
	
	dv_new(T0);
	for(i=0;i<5;i++)dv_new(t[i]);
	
	for(i=nb-1;i>=0;i--){
		int b=bn_get_bit(m,i);	
		for (int j = 0; j <= 5; j++)//5S,5M
		{
			fp_sqr(S1[j],V1[j + 1]);
			fp_mul(P1[j],V1[j],V1[j + 2]);
		}
		
		if (b == 0)
		{
			fp_add(k[0],P1[0],P1[1]);
			fp_sub(k[5],S1[0],S1[1]);
			fp_muln_low(t[0],k[0],k[5]);

			fp_sub(k[1],P1[0],P1[1]);
			fp_add(k[5],S1[0],S1[1]);
			fp_muln_low(t[1],k[1],k[5]);

			fp_add(k[2],P1[0],P1[2]);
			fp_sub(k[5],S1[0],S1[2]);
			fp_muln_low(t[2],k[2],k[5]);

			fp_sub(k[3],P1[0],P1[2]);
			fp_add(k[5],S1[0],S1[2]);
			fp_muln_low(t[3],k[3],k[5]);

			fp_sub(k[4],P1[1],P1[2]);
			fp_add(k[5],S1[1],S1[2]);
			fp_muln_low(t[4],k[4],k[5]);
			fp_addc_low(t[4],t[4],t[4]);

			fp_subc_low(T0,t[0],t[1]);
			fp_rdc(V1[0],T0);
			fp_subc_low(T0,t[2],t[3]);
			fp_rdc(V1[1],T0);
			fp_mul(V1[1],V1[1],alpha);


			fp_addc_low(T0,t[2],t[3]);
			fp_addc_low(t[2],t[0],t[1]);
			fp_addc_low(t[2],t[2],t[4]);
			fp_subc_low(T0,T0,t[2]);
			fp_rdc(V1[2],T0);

			fp_muln_low(t[0],S1[1],P1[3]);
			fp_muln_low(t[1],S1[3],P1[1]);
			fp_subc_low(T0,t[0],t[1]);
			fp_addc_low(T0,T0,T0);
			fp_rdc(V1[3],T0);
			fp_mul(V1[3],V1[3],alpha);


			fp_add(k[0],P1[2],P1[3]);
			fp_sub(k[5],S1[2],S1[3]);
			fp_muln_low(t[0],k[0],k[5]);

			fp_sub(k[1],P1[2],P1[3]);
			fp_add(k[5],S1[2],S1[3]);
			fp_muln_low(t[1],k[1],k[5]);

			fp_add(k[2],P1[2],P1[4]);
			fp_sub(k[5],S1[2],S1[4]);
			fp_muln_low(t[2],k[2],k[5]);

			fp_sub(k[3],P1[2],P1[4]);
			fp_add(k[5],S1[2],S1[4]);
			fp_muln_low(t[3],k[3],k[5]);

			fp_sub(k[4],P1[3],P1[4]);
			fp_add(k[5],S1[3],S1[4]);
			fp_muln_low(t[4],k[4],k[5]);
			fp_addc_low(t[4],t[4],t[4]);

			fp_subc_low(T0,t[0],t[1]);
			fp_rdc(V1[4],T0);
			fp_subc_low(T0,t[2],t[3]);
			fp_rdc(V1[5],T0);
			fp_mul(V1[5],V1[5],alpha);
			fp_addc_low(T0,t[2],t[3]);
			fp_addc_low(t[2],t[0],t[1]);
			fp_addc_low(t[2],t[2],t[4]);
			fp_subc_low(T0,T0,t[2]);
			fp_rdc(V1[6],T0);
			fp_muln_low(t[0],S1[3],P1[5]);
			fp_muln_low(t[1],S1[5],P1[3]);
			fp_subc_low(T0,t[0],t[1]);
			fp_addc_low(T0,T0,T0);
			fp_rdc(V1[7],T0);
			fp_mul(V1[7],V1[7],alpha);

		}
		else
		{
			fp_muln_low(t[0],S1[0],P1[2]);
			fp_muln_low(t[1],S1[2],P1[0]);
			fp_subc_low(T0,t[0],t[1]);
			fp_addc_low(T0,T0,T0);
			fp_rdc(V1[0],T0);
			fp_mul(V1[0],V1[0],alpha);

			fp_add(k[0],P1[1],P1[2]);
			fp_sub(k[5],S1[1],S1[2]);
			fp_muln_low(t[0],k[0],k[5]);

			fp_sub(k[1],P1[1],P1[2]);
			fp_add(k[5],S1[1],S1[2]);
			fp_muln_low(t[1],k[1],k[5]);

			fp_add(k[2],P1[1],P1[3]);
			fp_sub(k[5],S1[1],S1[3]);
			fp_muln_low(t[2],k[2],k[5]);

			fp_sub(k[3],P1[1],P1[3]);
			fp_add(k[5],S1[1],S1[3]);
			fp_muln_low(t[3],k[3],k[5]);

			fp_sub(k[4],P1[2],P1[3]);
			fp_add(k[5],S1[2],S1[3]);
			fp_muln_low(t[4],k[4],k[5]);
			fp_addc_low(t[4],t[4],t[4]);

			fp_subc_low(T0,t[0],t[1]);
			fp_rdc(V1[1],T0);
			fp_subc_low(T0,t[2],t[3]);
			fp_rdc(V1[2],T0);
			fp_mul(V1[2],V1[2],alpha);
			fp_addc_low(T0,t[2],t[3]);
			fp_addc_low(t[2],t[0],t[1]);
			fp_addc_low(t[2],t[2],t[4]);
			fp_subc_low(T0,T0,t[2]);
			fp_rdc(V1[3],T0);
			fp_muln_low(t[0],S1[2],P1[4]);
			fp_muln_low(t[1],S1[4],P1[2]);
			fp_subc_low(T0,t[0],t[1]);
			fp_addc_low(T0,T0,T0);
			fp_rdc(V1[4],T0);
			fp_mul(V1[4],V1[4],alpha);

			fp_add(k[0],P1[3],P1[4]);
			fp_sub(k[5],S1[3],S1[4]);
			fp_muln_low(t[0],k[0],k[5]);

			fp_sub(k[1],P1[3],P1[4]);
			fp_add(k[5],S1[3],S1[4]);
			fp_muln_low(t[1],k[1],k[5]);

			fp_add(k[2],P1[3],P1[5]);
			fp_sub(k[5],S1[3],S1[5]);
			fp_muln_low(t[2],k[2],k[5]);

			fp_sub(k[3],P1[3],P1[5]);
			fp_add(k[5],S1[3],S1[5]);
			fp_muln_low(t[3],k[3],k[5]);

			fp_sub(k[4],P1[4],P1[5]);
			fp_add(k[5],S1[4],S1[5]);
			fp_muln_low(t[4],k[4],k[5]);
			fp_addc_low(t[4],t[4],t[4]);

			fp_subc_low(T0,t[0],t[1]);
			fp_rdc(V1[5],T0);
			fp_subc_low(T0,t[2],t[3]);
			fp_rdc(V1[6],T0);
			fp_mul(V1[6],V1[6],alpha);
			fp_addc_low(T0,t[2],t[3]);
			fp_addc_low(t[2],t[0],t[1]);
			fp_addc_low(t[2],t[2],t[4]);
			fp_subc_low(T0,T0,t[2]);
			fp_rdc(V1[7],T0);
		}
	}
	
	dv_free(T0);
	for(i=0;i<5;i++)dv_free(t[i]);
}

void scalmul_net(ep_t mP,ep_t P,fp_t A,fp_t B,bn_t m,int version){
	int i;
	fp_t x_sp,y_sp,xp,yp,V1[8],alpha,S1[6],P1[6],k[6];
	fp_t t0,t1,t2;

/*========================================================*/
/*                        Init                            */
/*========================================================*/

	fp_null(x_sp);
	fp_null(y_sp);
	fp_null(xp);
	fp_null(yp);
	fp_null(alpha);
	fp_null(t0);
	fp_null(t1);
	fp_null(t2);
	for(i=0;i<8;i++){
		fp_null(V1[i]);
	}
	for(i=0;i<6;i++){
		fp_null(S1[i]);
		fp_null(P1[i]);
		fp_null(k[i]);
	}	
/*========================================================*/
	fp_new(x_sp);
	fp_new(y_sp);
	fp_new(xp);
	fp_new(yp);
	fp_new(alpha);
	fp_new(t0);
	fp_new(t1);
	fp_new(t2);
	for(i=0;i<8;i++){
		fp_new(V1[i]);
	}
	for(i=0;i<6;i++){
		fp_new(S1[i]);
		fp_new(P1[i]);
		fp_new(k[i]);
	}
/*========================================================*/
	fp_copy(xp,P->x);//xp
	fp_copy(yp,P->y);//yp

/*========================================================*/
/*                      Init data                         */
/*========================================================*/

/* fill out the first vector of the block */
	fp_set_dig(V1[3],1);//V1[3]=1
	fp_dbl(V1[4],yp);//the origin w(2,0)=2*yp
	/* the origin w(3,0)=3*xp^4+6*A*xp^2+12*B*xp-A*A */
	fp_sqr(t0,xp);//t0=xp^2
	fp_dbl(V1[5],t0);//2*t0
	fp_add(V1[5],V1[5],t0);//3*t0
	fp_dbl(t1,A);//2*A
	fp_add(t2,t1,t0);//2*A+xp^2
	fp_mul(V1[5],V1[5],t2);
	fp_sqr(t2,A);//A*A
	fp_sub(V1[5],V1[5],t2);
	fp_mul_dig(V1[6],B,12);
	fp_mul(V1[6],V1[6],xp);
	fp_add(V1[5],V1[5],V1[6]);
	/*the origin w(4,0) = 4*yp*(xp^6+5A*xp^4+20Bxp^3-5A^2xp^2-4AB*xp-8B^2-A^3)*/
	fp_add(t1,t1,t1);//4A
	fp_add(t1,t1,A);//5*A
	fp_mul(V1[6],t0,t1);//5Axp^2
	fp_sub(V1[2],t0,A);//xp^2-A
	fp_mul(V1[6],V1[6],V1[2]);
	fp_mul(t1,t0,xp);//xp^3
	fp_mul_dig(V1[2],B,20);
	fp_mul(V1[2],t1,V1[2]);//20B*xp^3
	fp_add(V1[6],V1[6],V1[2]);
	fp_sqr(t1,t1);//xp^6
	fp_add(V1[6],V1[6],t1);
	fp_mul_dig(t1,A,4);
	fp_mul(t1,t1,B);
	fp_mul(t1,t1,xp);//4AB*xp
	fp_sub(V1[6],V1[6],t1);
	fp_sqr(t1,B);
	fp_mul_dig(t1,t1,8);
	fp_sub(V1[6],V1[6],t1);
	fp_mul(t1,t2,A);
	fp_sub(V1[6],V1[6],t1);
	//8*B*B=0
	fp_dbl(t1,yp);//2*yp
	fp_add(t2,t1,t1);//4*yp
	fp_mul(V1[6],V1[6],t2);
	/*the orginal W(5,0)*/
	fp_sqr(V1[7],V1[5]);
	fp_mul(V1[7],V1[7],V1[5]);
	fp_sqr(t1,V1[4]);
	fp_mul(t1,t1,V1[4]);
	fp_mul(t1,t1,V1[6]);
	fp_sub(V1[7],t1,V1[7]);

	fp_zero(V1[2]);//the orginal W(0,0)=0
	fp_neg(V1[1],V1[3]);//the orginal W(-1,0)
	fp_neg(V1[0],V1[4]);//the orginal W(-2,0)
	
	fp_dbl(t0,yp);//2*yp
	fp_inv(alpha,t0);//1/2*yp
	switch(version){
		case 1:
			DoubleAndAdd01(V1,S1,P1,k,alpha,m);
			break;
		case 2:
			DoubleAndAdd01_lz(V1,S1,P1,k,alpha,m);
			break;
		case 3:
			DoubleAndAdd02(V1,S1,P1,k,alpha,m);
			break;
		case 4:
			DoubleAndAdd02_lz(V1,S1,P1,k,alpha,m);
			break;
		default:
			DoubleAndAdd01(V1,S1,P1,k,alpha,m);
			break;
	}

/*x_sp = xp - (V1[2] * V1[4]) / (V1[3] * V1[3])*/
	fp_sqr(y_sp,V1[3]);
	fp_inv(t0,y_sp);
	fp_mul(t1,V1[2],V1[4]);	
	fp_mul(t1,t1,t0);
	fp_sub(x_sp,xp,t1);//notice that when m=r,V1[3]=0
	fp_copy(mP->x,x_sp);
/*y_sp=(V1[2]^2*V1[5]-V1[4]^2*V1[1])/4*yq*V1[3]^3*/
	fp_mul(t0,y_sp,V1[3]);
	fp_mul(y_sp,t2,t0);
	fp_inv(y_sp,y_sp);
	fp_sqr(t0,V1[2]);
	fp_mul(t0,t0,V1[5]);
	fp_sqr(t1,V1[4]);
	fp_mul(t1,t1,V1[1]);
	fp_sub(t0,t0,t1);
	fp_mul(y_sp,t0,y_sp);
	fp_copy(mP->y,y_sp);	
	fp_set_dig(mP->z,1);

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	fp_free(x_sp);
	fp_free(y_sp);
	fp_free(xp);
	fp_free(yp);
	fp_free(alpha);
	fp_free(t0);
	fp_free(t1);
	fp_free(t2);
	for(i=0;i<8;i++){
		fp_free(V1[i]);
	}
	for(i=0;i<6;i++){
		fp_free(S1[i]);
		fp_free(P1[i]);
		fp_free(k[i]);
	}
}
