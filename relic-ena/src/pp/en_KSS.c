/*==========================================================================*/
/*===========                                              =================*/
/*               Original Elliptic Net Alg on KSS curves                    */
/*===========                                              =================*/
/*==========================================================================*/

#include"ena.h"

static void precomp(fp18_t F,fp18_t x_q1,fp18_t y_q1,fp18_t xq,fp18_t yq,fp18_t initial_data[],fp18_t V2[],fp18_t V1[]){
	fp18_t alpha,beta,gamma1,S_0,P_0,t0;
	fp18_t S1[6],P1[6],tv2[3],tv1[8];
	int i;

/*========================================================*/
/*                          Init                          */
/*========================================================*/
	TRY{
		fp18_null(alpha);
		fp18_null(beta);
		fp18_null(gamma1);
		fp18_null(S_0);
		fp18_null(P_0);
		fp18_null(t0);

		for(i=0;i<6;i++){
			fp18_null(S1[i]);
			fp18_null(P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_null(tv2[i]);
		}
		for(i=0;i<8;i++){
			fp18_null(tv1[i]);
		}
/*========================================================*/	
		fp18_new(alpha);
		fp18_new(beta);
		fp18_new(gamma1);
		fp18_new(S_0);
		fp18_new(P_0);
		fp18_new(t0);

		for(i=0;i<6;i++){
			fp18_new(S1[i]);
			fp18_new(P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_new(tv2[i]);
		}
		for(i=0;i<8;i++){
			fp18_new(tv1[i]);
		}

/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

		for(i=0;i<3;i++){
			fp18_copy(tv2[i],V2[i]);
		}
		for(i=0;i<8;i++){
			fp18_copy(tv1[i],V1[i]);
		}
		
		
		//initial data contains the precomputed inverses
		fp18_copy(alpha,initial_data[0]);		
		fp18_copy(beta,initial_data[1]);// inverse of W(-1,1)
		fp18_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
		fp18_sqr(S_0,tv2[1]);//1S_12
		fp18_mul(P_0,tv2[0],tv2[2]);//1M_12
		for (int j = 0; j <= 5; j++)//6S_12,6M_12
		{
			fp18_sqr(S1[j],tv1[j + 1]);
			fp18_mul(P1[j],tv1[j],tv1[j + 2]);
		}
		for (int j = 1; j <= 4; j++)
		{
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
			fp18_mul(tv1[2*j-2],S1[j-1],P1[j+1]);
			fp18_mul(t0,S1[j+1],P1[j-1]);
			fp18_sub(tv1[2*j-2],tv1[2*j-2],t0);
			fp18_mul_dxs(tv1[2*j-2],tv1[2*j-2],alpha);
			/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
			fp18_mul(tv1[2*j-1],S1[j],P1[j+1]);		
			fp18_mul(t0,S1[j+1],P1[j]);
			fp18_sub(tv1[2*j-1],tv1[2*j-1],t0);
		}

		/* tv2[0] = S1[2] * P_0 - P1[2] * S_0 */
		fp18_mul(tv2[0],S1[2],P_0);
		fp18_mul(t0,P1[2],S_0);
		fp18_sub(tv2[0],tv2[0],t0);
		/* tv2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
		fp18_mul(tv2[1],S1[3],P_0);
		fp18_mul(t0,P1[3],S_0);
		fp18_sub(tv2[1],tv2[1],t0);
		fp18_mul_dxs(tv2[1],tv2[1],beta);
		/* tv2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
		fp18_mul(tv2[2],P1[4],S_0);
		fp18_mul(t0,S1[4],P_0);
		fp18_sub(tv2[2],tv2[2],t0);
		fp18_mul(tv2[2],tv2[2],gamma1);

	/*x_q1 = xq - (tv1[2] * tv1[4]) / (tv1[3] * tv1[3])*/
		fp18_sqr(t0,tv1[3]);
		fp18_inv(x_q1,t0);
		fp18_mul(y_q1,tv1[2],tv1[4]);	
		fp18_mul(y_q1,y_q1,x_q1);
		fp18_sub(x_q1,xq,y_q1);//notice that when m=r,tv1[3]=0
		
	/*y_q1=(tv1[2]^2*tv1[5]-tv1[4]^2*tv1[1])/4*yq*tv1[3]^3*/
		fp18_mul(y_q1,t0,tv1[3]);
		fp6_zero(t0[1]);
		fp6_zero(t0[2]);
		fp6_dbl(t0[0],yq[0]);
		fp6_dbl(t0[0],t0[0]);
		fp18_mul(y_q1,y_q1,t0);
		fp18_inv(y_q1,y_q1);
		fp18_sqr(t0,tv1[2]);
		fp18_mul(t0,t0,tv1[5]);
		fp18_sqr(gamma1,tv1[4]);
		fp18_mul(gamma1,gamma1,tv1[1]);
		fp18_sub(t0,t0,gamma1);
		fp18_mul(y_q1,t0,y_q1);

		fp18_copy(F,tv2[1]);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}	

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {

		fp18_free(alpha);	
		fp18_free(beta);
		fp18_free(gamma1);
		fp18_free(S_0);
		fp18_free(P_0);
		fp18_free(t0);

		for(i=0;i<6;i++){
			fp18_free(S1[i]);
			fp18_free(P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_free(tv2[i]);
		}
		for(i=0;i<8;i++){
			fp18_free(tv1[i]);
		}
	}
}


static void precomp_lz(fp18_t F,fp18_t x_q1,fp18_t y_q1,fp18_t xq,fp18_t yq,fp18_t initial_data[],fp18_t V2[],fp18_t V1[]){
	fp18_t alpha,beta,gamma1,S_0,P_0;
	fp18_t S1[6],P1[6],tv2[3],tv1[8];
	dv18_t t0,t1;
	int i;

/*========================================================*/
/*                          Init                          */
/*========================================================*/
	TRY{
		fp18_null(alpha);
		fp18_null(beta);
		fp18_null(gamma1);
		fp18_null(S_0);
		fp18_null(P_0);
		dv18_null(t0);
		dv18_null(t1);

		for(i=0;i<6;i++){
			fp18_null(S1[i]);
			fp18_null(P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_null(tv2[i]);
		}
		for(i=0;i<8;i++){
			fp18_null(tv1[i]);
		}
/*========================================================*/	
		fp18_new(alpha);
		fp18_new(beta);
		fp18_new(gamma1);
		fp18_new(S_0);
		fp18_new(P_0);
		dv18_new(t0);
		dv18_new(t1);

		for(i=0;i<6;i++){
			fp18_new(S1[i]);
			fp18_new(P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_new(tv2[i]);
		}
		for(i=0;i<8;i++){
			fp18_new(tv1[i]);
		}

/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

		for(i=0;i<3;i++){
			fp18_copy(tv2[i],V2[i]);
		}
		for(i=0;i<8;i++){
			fp18_copy(tv1[i],V1[i]);
		}
		
		
		//initial data contains the precomputed inverses
		fp18_copy(alpha,initial_data[0]);		
		fp18_copy(beta,initial_data[1]);// inverse of W(-1,1)
		fp18_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
		fp18_sqr(S_0,tv2[1]);//1S_12
		fp18_mul(P_0,tv2[0],tv2[2]);//1M_12
		for (int j = 0; j <= 5; j++)//6S_12,6M_12
		{
			fp18_sqr(S1[j],tv1[j + 1]);
			fp18_mul(P1[j],tv1[j],tv1[j + 2]);
		}
		for (int j = 1; j <= 4; j++)
		{
		/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
			fp18_muln_low(t0,S1[j-1],P1[j+1]);
			fp18_muln_low(t1,S1[j+1],P1[j-1]);
			fp18_subc_low(t0,t0,t1);
			fp18_rdc(tv1[2*j-2],t0);
			fp18_mul_dxs(tv1[2*j-2],tv1[2*j-2],alpha);
		/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
			fp18_muln_low(t0,S1[j],P1[j+1]);		
			fp18_muln_low(t1,S1[j+1],P1[j]);
			fp18_subc_low(t0,t0,t1);
			fp18_rdc(tv1[2*j-1],t0);
		}

		/* tv2[0] = S1[2] * P_0 - P1[2] * S_0 */
		fp18_muln_low(t0,S1[2],P_0);
		fp18_muln_low(t1,P1[2],S_0);
		fp18_subc_low(t0,t0,t1);
		fp18_rdc(tv2[0],t0);
		/* tv2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
		fp18_muln_low(t0,S1[3],P_0);
		fp18_muln_low(t1,P1[3],S_0);
		fp18_subc_low(t0,t0,t1);
		fp18_rdc(tv2[1],t0);
		fp18_mul_dxs(tv2[1],tv2[1],beta);
		/* tv2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
		fp18_muln_low(t0,P1[4],S_0);
		fp18_muln_low(t1,S1[4],P_0);
		fp18_subc_low(t0,t0,t1);
		fp18_rdc(tv2[2],t0);
		fp18_mul(tv2[2],tv2[2],gamma1);

	/*x_q1 = xq - (tv1[2] * tv1[4]) / (tv1[3] * tv1[3])*/
		fp18_sqr(beta,tv1[3]);
		fp18_inv(x_q1,beta);
		fp18_mul(y_q1,tv1[2],tv1[4]);	
		fp18_mul(y_q1,y_q1,x_q1);
		fp18_sub(x_q1,xq,y_q1);//notice that when m=r,tv1[3]=0
		
	/*y_q1=(tv1[2]^2*tv1[5]-tv1[4]^2*tv1[1])/4*yq*tv1[3]^3*/
		fp18_mul(y_q1,beta,tv1[3]);
		fp6_zero(beta[1]);
		fp6_zero(beta[2]);
		fp6_dbl(beta[0],yq[0]);
		fp6_dbl(beta[0],beta[0]);
		fp18_mul(y_q1,y_q1,beta);
		fp18_inv(y_q1,y_q1);
		fp18_sqr(beta,tv1[2]);
		fp18_mul(beta,beta,tv1[5]);
		fp18_sqr(gamma1,tv1[4]);
		fp18_mul(gamma1,gamma1,tv1[1]);
		fp18_sub(beta,beta,gamma1);
		fp18_mul(y_q1,beta,y_q1);

		fp18_copy(F,tv2[1]);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}	

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {

		fp18_free(alpha);	
		fp18_free(beta);
		fp18_free(gamma1);
		fp18_free(S_0);
		fp18_free(P_0);
		dv18_free(t0);
		dv18_free(t1);

		for(i=0;i<6;i++){
			fp18_free(S1[i]);
			fp18_free(P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_free(tv2[i]);
		}
		for(i=0;i<8;i++){
			fp18_free(tv1[i]);
		}
	}
}
/*origin double and add alg*/
static void DoubleAndAdd(fp18_t initial_data[],fp18_t V2[],fp18_t V1[],bn_t m){
	fp18_t alpha,beta,gamma1,delta,S_0,P_0,t0;
	fp18_t S1[6],P1[6];
	int nb,i;

/*========================================================*/
/*                          Init                          */
/*========================================================*/
TRY{
	fp18_null(alpha);
	fp18_null(beta);
	fp18_null(gamma1);
	fp18_null(delta);
	fp18_null(S_0);
	fp18_null(P_0);
	fp18_null(t0);

	for(i=0;i<6;i++){
		fp18_null(S1[i]);
		fp18_null(P1[i]);
	}
/*========================================================*/	
	fp18_new(alpha);
	fp18_new(beta);
	fp18_new(gamma1);
	fp18_new(delta);
	fp18_new(S_0);
	fp18_new(P_0);
	fp18_new(t0);

	for(i=0;i<6;i++){
		fp18_new(S1[i]);
		fp18_new(P1[i]);
	}

/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

	nb = bn_bits(m);
	nb = nb - 1;
	
	//initial data contains the precomputed inverses
	fp18_copy(alpha,initial_data[0]);		
	fp18_copy(beta,initial_data[1]);// inverse of W(-1,1)
	fp18_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
	fp18_copy(delta,initial_data[4]); // inverse of W(1,1)	
	for(i=nb-1;i>=0;i--){
		int b=bn_get_bit(m,i);	
		fp18_sqr(S_0,V2[1]);//1S_12
		fp18_mul(P_0,V2[0],V2[2]);//1M_12
		for (int j = 0; j <= 5; j++)//6S_12,6M_12
		{
			fp18_sqr(S1[j],V1[j + 1]);
			fp18_mul(P1[j],V1[j],V1[j + 2]);
		}
		
		if (b == 0)
		{
			for (int j = 1; j <= 4; j++)//j=3 or 4
			{
			/* S[j - 1] * P1[j] - S[j] * P1[j - 1] */
				fp18_mul(V1[2*j-2],S1[j-1],P1[j]);
				fp18_mul(t0,S1[j],P1[j-1]);
				fp18_sub(V1[2*j-2],V1[2*j-2],t0);
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
				fp18_mul(V1[2*j-1],S1[j-1],P1[j+1]);	
				fp18_mul(t0,S1[j+1],P1[j-1]);
				fp18_sub(V1[2*j-1],V1[2*j-1],t0);
				fp18_mul_dxs(V1[2*j-1],V1[2*j-1],alpha);
			}

			/* V2[0] = (S1[1] * P_0 - P1[1] * S_0) * delta*/
			fp18_mul(V2[0],S1[1],P_0);
			fp18_mul(t0,P1[1],S_0);
			fp18_sub(V2[0],V2[0],t0);
			fp18_mul(V2[0],V2[0],delta);
			/*V2[1] = S1[2] * P_0 - P1[2] * S_0*/
			fp18_mul(V2[1],S1[2],P_0);
			fp18_mul(t0,P1[2],S_0);
			fp18_sub(V2[1],V2[1],t0);
			/* V2[2] = (S1[3] * P_0 - P1[3] * S_0) * beta */
			fp18_mul(V2[2],S1[3],P_0);
			fp18_mul(t0,P1[3],S_0);
			fp18_sub(V2[2],V2[2],t0);
			fp18_mul_dxs(V2[2],V2[2],beta);
		}
		else
		{
			for (int j = 1; j <= 4; j++)
			{
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
				fp18_mul(V1[2*j-2],S1[j-1],P1[j+1]);
				fp18_mul(t0,S1[j+1],P1[j-1]);
				fp18_sub(V1[2*j-2],V1[2*j-2],t0);
				fp18_mul_dxs(V1[2*j-2],V1[2*j-2],alpha);
			/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
				fp18_mul(V1[2*j-1],S1[j],P1[j+1]);		
				fp18_mul(t0,S1[j+1],P1[j]);
				fp18_sub(V1[2*j-1],V1[2*j-1],t0);
			}

			/* V2[0] = S1[2] * P_0 - P1[2] * S_0 */
			fp18_mul(V2[0],S1[2],P_0);
			fp18_mul(t0,P1[2],S_0);
			fp18_sub(V2[0],V2[0],t0);
			/* V2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
			fp18_mul(V2[1],S1[3],P_0);
			fp18_mul(t0,P1[3],S_0);
			fp18_sub(V2[1],V2[1],t0);
			fp18_mul_dxs(V2[1],V2[1],beta);
			/* V2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
			fp18_mul(V2[2],P1[4],S_0);
			fp18_mul(t0,S1[4],P_0);
			fp18_sub(V2[2],V2[2],t0);
			fp18_mul(V2[2],V2[2],gamma1);
		}
	}

}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}	

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {

	fp18_free(alpha);	
	fp18_free(beta);
	fp18_free(gamma1);
	fp18_free(delta);
	fp18_free(S_0);
	fp18_free(P_0);
	fp18_free(t0);

	for(i=0;i<6;i++){
		fp18_free(S1[i]);
		fp18_free(P1[i]);
	}
}
}

/*origin double and add alg with lazy reduction*/
static void DoubleAndAdd_lz(fp18_t initial_data[],fp18_t V2[],fp18_t V1[],bn_t m){
	fp18_t alpha,beta,gamma1,delta,S_0,P_0;
	fp18_t S1[6],P1[6];
	dv18_t t0,t1;
	int nb,i;

/*========================================================*/
/*                        Init                            */
/*========================================================*/
TRY{
	fp18_null(alpha);
	fp18_null(beta);
	fp18_null(gamma1);
	fp18_null(delta);
	fp18_null(S_0);
	fp18_null(P_0);
	dv18_null(t0);
	dv18_null(t1);

	for(i=0;i<6;i++){
		fp18_null(S1[i]);
		fp18_null(P1[i]);
	}
/*========================================================*/
	fp18_new(alpha);
	fp18_new(beta);
	fp18_new(gamma1);
	fp18_new(delta);
	fp18_new(S_0);
	fp18_new(P_0);
	dv18_new(t0);
	dv18_new(t1);

	for(i=0;i<6;i++){
		fp18_new(S1[i]);
		fp18_new(P1[i]);
	}

/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

	nb = bn_bits(m);
	nb = nb - 1;
	
	//initial data contains the precomputed inverses
	fp18_copy(alpha,initial_data[0]);		
	fp18_copy(beta,initial_data[1]);// inverse of W(-1,1)
	fp18_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
	fp18_copy(delta,initial_data[4]); // inverse of W(1,1)	
	for(i=nb-1;i>=0;i--){
		int b=bn_get_bit(m,i);	
		fp18_sqr(S_0,V2[1]);//1S_12
		fp18_mul(P_0,V2[0],V2[2]);//1M_12
		for (int j = 0; j <= 5; j++)//6S_12,6M_12
		{
			fp18_sqr(S1[j],V1[j + 1]);
			fp18_mul(P1[j],V1[j],V1[j + 2]);
		}
		
		if (b == 0)
		{
			for (int j = 1; j <= 4; j++)//j=3 or 4
			{
			/* S[j - 1] * P1[j] - S[j] * P1[j - 1] */
				fp18_muln_low(t0,S1[j-1],P1[j]);
				fp18_muln_low(t1,S1[j],P1[j-1]);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V1[2*j-2],t0);
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
				fp18_muln_low(t0,S1[j-1],P1[j+1]);	
				fp18_muln_low(t1,S1[j+1],P1[j-1]);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V1[2*j-1],t0);
				fp18_mul_dxs(V1[2*j-1],V1[2*j-1],alpha);
			}

			/* V2[0] = (S1[1] * P_0 - P1[1] * S_0) * delta*/
			fp18_muln_low(t0,S1[1],P_0);
			fp18_muln_low(t1,P1[1],S_0);
			fp18_subc_low(t0,t0,t1);
			fp18_rdc(V2[0],t0);
			fp18_mul(V2[0],V2[0],delta);
			/*V2[1] = S1[2] * P_0 - P1[2] * S_0*/
			fp18_muln_low(t0,S1[2],P_0);
			fp18_muln_low(t1,P1[2],S_0);
			fp18_subc_low(t0,t0,t1);
			fp18_rdc(V2[1],t0);
			/* V2[2] = (S1[3] * P_0 - P1[3] * S_0) * beta */
			fp18_muln_low(t0,S1[3],P_0);
			fp18_muln_low(t1,P1[3],S_0);
			fp18_subc_low(t0,t0,t1);
			fp18_rdc(V2[2],t0);
			fp18_mul_dxs(V2[2],V2[2],beta);
		}
		else
		{
			for (int j = 1; j <= 4; j++)
			{
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
				fp18_muln_low(t0,S1[j-1],P1[j+1]);
				fp18_muln_low(t1,S1[j+1],P1[j-1]);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V1[2*j-2],t0);
				fp18_mul_dxs(V1[2*j-2],V1[2*j-2],alpha);
			/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
				fp18_muln_low(t0,S1[j],P1[j+1]);		
				fp18_muln_low(t1,S1[j+1],P1[j]);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V1[2*j-1],t0);
			}

			/* V2[0] = S1[2] * P_0 - P1[2] * S_0 */
			fp18_muln_low(t0,S1[2],P_0);
			fp18_muln_low(t1,P1[2],S_0);
			fp18_subc_low(t0,t0,t1);
			fp18_rdc(V2[0],t0);
			/* V2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
			fp18_muln_low(t0,S1[3],P_0);
			fp18_muln_low(t1,P1[3],S_0);
			fp18_subc_low(t0,t0,t1);
			fp18_rdc(V2[1],t0);
			fp18_mul_dxs(V2[1],V2[1],beta);
			/* V2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
			fp18_muln_low(t0,P1[4],S_0);
			fp18_muln_low(t1,S1[4],P_0);
			fp18_subc_low(t0,t0,t1);
			fp18_rdc(V2[2],t0);
			fp18_mul(V2[2],V2[2],gamma1);
		}
	}

}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}	

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {

	fp18_free(alpha);	
	fp18_free(beta);
	fp18_free(gamma1);
	fp18_free(delta);
	fp18_free(S_0);
	fp18_free(P_0);
	dv18_free(t0);
	dv18_free(t1);	

	for(i=0;i<6;i++){
		fp18_free(S1[i]);
		fp18_free(P1[i]);
	}
}
	
}

void opt_en_kss(fp18_t F,ep3_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version){
	int i=0,sign=0;
	fp18_t initial_data[7],V2[3],V1[8],xq,yq,w_n1,w_2n,w_21,w_11,k0,k1,B0,x_q1,y_q1,x_mq,y_mq;
	fp_t xp,yp;
	fp3_t a,d;	
	fp6_t _t0,_t1,_t2;
	bn_t l;

/*========================================================*/
/*                        Init                            */
/*========================================================*/
	TRY{
		bn_null(l);
		fp18_null(k0);
		fp18_null(k1);	
		fp18_null(xq);
		fp18_null(yq);
		fp18_null(w_n1);
		fp18_null(w_2n);
		fp18_null(w_21);
		fp18_null(w_11);
		fp18_null(B0);
		for(i=0;i<7;i++)	
			fp18_null(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_null(V2[i]);
		}
		fp3_null(a)
		fp3_null(d);
		fp_null(xp);
		fp_null(yp);
		fp18_null(x_q1);
		fp18_null(y_q1);
		fp18_null(x_mq);
		fp18_null(y_mq);
		fp6_null(_t0);
		fp6_null(_t1);
		fp6_null(_t2);
		for(i=0;i<8;i++){
			fp18_null(V1[i]);
		}	
/*========================================================*/	
		bn_new(l);
		fp18_new(k0);
		fp18_new(k1);	
		fp18_new(xq);
		fp18_new(yq);
		fp18_new(w_n1);
		fp18_new(w_2n);
		fp18_new(w_21);
		fp18_new(w_11);
		fp18_new(B0);
		for(i=0;i<7;i++)	
			fp18_new(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_new(V2[i]);
		}
		fp3_new(a);
		fp3_new(d);
		fp_new(xp);
		fp_new(yp);
		fp18_new(x_q1);
		fp18_new(y_q1);
		fp18_new(x_mq);
		fp18_new(y_mq);
		fp6_new(_t0);
		fp6_new(_t1);
		fp6_new(_t2);
		for(i=0;i<8;i++){
			fp18_new(V1[i]);
		}
/*========================================================*/
		fp18_set_dig(F,1);
		if (ep_is_infty(P) || ep2_is_infty(Q)) return;
		fp3_copy(a,Q->x);//a=xq
		fp3_copy(d,Q->y);//d=yq
		fp_copy(xp,P->x);//xp
		fp_copy(yp,P->y);//yp

		fp18_zero(xq);
		fp18_zero(yq);
		fp_copy(xq[2][0][0],a[0]);
		fp_copy(xq[2][2][0],a[1]);
		fp_copy(xq[2][1][1],a[2]);
		fp_copy(yq[0][1][0],d[0]);
		fp_copy(yq[0][0][1],d[1]);
		fp_copy(yq[0][2][1],d[2]);	

		fp18_zero(B0);
		fp_copy(B0[0][0][0],B);
		if(bn_sign(m)==1){
			bn_neg(l,m);
			sign=1;	
		}
		else{
			bn_copy(l,m);
		}

/*========================================================*/
/*                      Init data                         */
/*========================================================*/

	/* fill out the first vector of the block */

		fp18_set_dig(V1[3],1);//V1[3]=1
		fp6_dbl(V1[4][0],yq[0]);//the origin w(2,0)=2*yq
		fp6_zero(V1[4][1]);
		fp6_zero(V1[4][2]);
		/* the origin w(3,0)=3*xq^4+12*B*xq=3*xq*(xq^3+4*B) */
		fp6_zero(V1[5][0]);
		fp6_zero(V1[5][1]);		
		fp6_dbl(V1[5][2],xq[2]);
		fp6_add(V1[5][2],V1[5][2],xq[2]);//3*xq
		fp18_sqr(k0,xq);//t0=xq^2
		fp18_mul(k0,k0,xq);//xq^3

		fp6_dbl(_t2,B0[0]);
		fp6_add(_t2,_t2,_t2);//4B
		fp6_zero(k1[1]);
		fp6_zero(k1[2]);
		fp6_copy(k1[0],_t2);//k1=4B
		fp18_add(V1[6],k0,k1);
		fp18_mul(V1[5],V1[5],V1[6]);
		/*the origin w(4,0)=4*yq*(xq^6+4B(5*xq^3-2*B))*/
		fp6_zero(V1[6][1]);	
		fp6_dbl(_t1,k0[0]);//2*xq^3
		fp6_add(_t1,_t1,_t1);
		fp6_add(_t1,_t1,k0[0]);//5*xq^3
		fp6_dbl(V1[6][0],B0[0]);//2B
		fp6_sub(V1[6][0],_t1,V1[6][0]);
		fp6_zero(V1[6][1]);
		fp6_zero(V1[6][2]);
		fp18_mul(V1[6],V1[6],k1);
		fp18_sqr(k1,k0);//xq^6
		fp18_add(V1[6],V1[6],k1);
		fp6_zero(V1[2][1]);
		fp6_zero(V1[2][2]);
		fp6_dbl(V1[2][0],V1[4][0]);
		fp6_copy(_t2,V1[2][0]);//t2=4*yq
		fp18_mul(V1[6],V1[6],V1[2]);
		/*the orginal W(5,0)*/
		fp18_sqr(V1[7],V1[5]);
		fp18_mul(V1[7],V1[5],V1[7]);
		fp18_sqr(V1[2],V1[4]);
		fp18_mul(V1[2],V1[2],V1[4]);
		fp18_mul(V1[2],V1[2],V1[6]);
		fp18_sub(V1[7],V1[2],V1[7]);

		fp18_zero(V1[2]);//the orginal W(0,0)=0
		fp18_neg(V1[1],V1[3]);//the orginal W(-1,0)
		fp18_neg(V1[0],V1[4]);//the orginal W(-2,0)

	/* fill out the second vector of the block */

		fp18_set_dig(V2[0],1);// the original W(0,1)
		fp18_set_dig(V2[1],1);// the original W(1,1)
		/* the original W(2,1) =2*xq + xp - ((yp - yq)/(xp-xq))^2*/
		fp18_zero(k0);
		fp6_neg(k0[2],xq[2]);
		fp_copy(k0[0][0][0],xp);
		fp18_inv(k0,k0);//1/(xp-xq)
		fp18_neg(k1,yq);
		fp_copy(k1[0][0][0],yp);
		fp18_mul(k1,k1,k0);
		fp18_sqr(k1,k1);

		fp6_dbl(k0[2],xq[2]);
		fp6_zero(k0[0]);
		fp6_zero(k0[1]);
		fp_copy(k0[0][0][0],xp);//k0=2*xq + xp
		fp18_sub(V2[2],k0,k1);

		fp18_copy(w_n1,xq);
		fp_neg(w_n1[0][0][0],xp);//w_n1=xq-xp
		fp18_sqr(k1,w_n1);//k1=w_n1^2
		fp18_mul(k0,k1,k0);
		fp18_neg(w_2n,k0);
		fp18_copy(k0,yq);
		fp_copy(k0[0][0][0],yp);
		fp18_sqr(k0,k0);//(yp+y_q2)^2
		fp18_add(w_2n,w_2n,k0);	
		fp18_copy(w_21,V2[2]);//the origin w(2,1)
		
		//V2[0]=V2[0]
		fp18_copy(V2[1],w_n1);//V2[1]=V2[1]*w_n1=1*w_n1=w_n1,new w(1,1)
		fp18_mul(V2[2],V2[2],k1);//new w(2,1)

		fp18_inv(k0,k1);//1/w_n1^2
		fp18_mul(w_2n,w_2n,k0);//w(-2,1)
		fp18_set_dig(w_n1,1);//w(1,-1)
		fp18_copy(w_11,V2[1]);//w(1,1)
		
		fp18_inv(initial_data[0],V1[4]);//w(2,0)^(-1)
		fp18_inv(initial_data[1],w_n1);//w_n1=1,w(-1,1)^(-1)
		fp18_inv(initial_data[2],w_2n);//w(-2,1)^(-1)
		fp18_inv(initial_data[3],w_21);//w(2,1)^(-1)
		fp18_inv(initial_data[4],w_11);//w(1,1)^(-1)
		fp18_sqr(initial_data[5],V1[4]);
		fp18_mul(initial_data[6],V1[3],V1[5]);
					
		switch(version){
			case 1:
				precomp(F,x_q1,y_q1,xq,yq,initial_data,V2,V1);
				DoubleAndAdd(initial_data,V2,V1,l);
				break;
			case 2:
				precomp_lz(F,x_q1,y_q1,xq,yq,initial_data,V2,V1);
				DoubleAndAdd_lz(initial_data,V2,V1,l);
				break;
			default:
				precomp(F,x_q1,y_q1,xq,yq,initial_data,V2,V1);
				DoubleAndAdd(initial_data,V2,V1,l);
				break;
		}
	//	fp18_inv(k0,V1[3]);
	//	fp18_mul(F,k0,V2[1]);
	/*x_mq = xq - (V1[2] * V1[4]) / (V1[3] * V1[3])*/
		fp18_sqr(k0,V1[3]);
		fp18_inv(x_mq,k0);
		fp18_mul(k1,V1[2],V1[4]);	
		fp18_mul(k1,k1,x_mq);
		fp18_sub(x_mq,xq,k1);//notice that when m=r,V1[3]=0
		
	/*y_mq=(V1[2]^2*V1[5]-V1[4]^2*V1[1])/4*yq*V1[3]^3*/
		fp18_mul(k0,k0,V1[3]);
		fp6_zero(k1[1]);
		fp6_zero(k1[2]);
		fp6_dbl(k1[0],yq[0]);
		fp6_dbl(k1[0],k1[0]);
		fp18_mul(y_mq,k1,k0);
		fp18_inv(y_mq,y_mq);
		fp18_sqr(k0,V1[2]);
		fp18_mul(k0,k0,V1[5]);
		fp18_sqr(k1,V1[4]);
		fp18_mul(k1,k1,V1[1]);
		fp18_sub(k0,k0,k1);
		fp18_mul(y_mq,k0,y_mq);

		/* f_{-a,Q}(P) = 1/f_{a,Q}(P). */
		if(sign==1){
			fp18_inv(k0, V2[1]);
			fp18_neg(y_mq,y_mq);
		}
		else{
			fp18_copy(k0,V2[1]);	
		}

		fp18_frb(F,F,1);//l^p
		fp18_frb(x_q1,x_q1,1);//x_q1^p
		fp18_frb(y_q1,y_q1,1);//y_q1^p
		fp18_mul(F,F,k0);
		
		fp18_zero(k0);
		fp_sub(k0[2][0][0],x_mq[2][0][0], x_q1[2][0][0]);
		fp_sub(k0[2][2][0],x_mq[2][2][0], x_q1[2][2][0]);
		fp_sub(k0[2][1][1],x_mq[2][1][1], x_q1[2][1][1]);	
		fp18_inv(k0,k0);
		fp18_copy(k1,y_q1);
		fp_sub(k1[0][1][0],y_mq[0][1][0], y_q1[0][1][0]);
		fp_sub(k1[0][0][1],y_mq[0][0][1], y_q1[0][0][1]);
		fp_sub(k1[0][2][1],y_mq[0][2][1], y_q1[0][2][1]);	
		fp18_mul_dxs(k0,k0,k1);
		fp_sub(x_q1[0][0][0],x_q1[0][0][0],xp);
		fp18_mul(k0,k0,x_q1);
		fp18_neg(k0,k0);
		fp_sub(y_q1[0][0][0],y_q1[0][0][0],yp);
		fp18_add(k0,k0,y_q1);
		fp18_mul(F,F,k0);
		
		pp_exp_k18(F,F);

	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}	

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {

		bn_free(l);
		fp18_free(k0);
		fp18_free(k1);	
		fp18_free(xq);
		fp18_free(yq);
		fp18_free(w_n1);
		fp18_free(w_2n);
		fp18_free(w_21);
		fp18_free(w_11);
		fp18_free(B0);
		for(i=0;i<7;i++)	
			fp18_free(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_free(V2[i]);
		}
		fp3_free(a);
		fp3_free(d);
		fp_free(xp);
		fp_free(yp);
		fp18_free(x_q1);
		fp18_free(y_q1);
		fp18_free(x_mq);
		fp18_free(y_mq);
		fp6_free(_t0);
		fp6_free(_t1);
		fp6_free(_t2);
		for(i=0;i<8;i++){
			fp18_free(V1[i]);
		}
	}

}

void opt_en_kss_lz(fp18_t F,ep3_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version){
	int i=0,sign=0;
	fp18_t initial_data[7],V2[3],V1[8],xq,yq,w_n1,w_2n,w_21,w_11,k0,k1,B0,x_q1,y_q1,x_mq,y_mq;
	fp_t xp,yp;
	fp3_t a,d;	
	dv18_t _t0,_t1;
	bn_t l;

/*========================================================*/
/*                        Init                            */
/*========================================================*/
	TRY{
		bn_null(l);
		fp18_null(k0);
		fp18_null(k1);	
		fp18_null(xq);
		fp18_null(yq);
		fp18_null(w_n1);
		fp18_null(w_2n);
		fp18_null(w_21);
		fp18_null(w_11);
		fp18_null(B0);
		for(i=0;i<7;i++)	
			fp18_null(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_null(V2[i]);
		}
		fp3_null(a)
		fp3_null(d);
		fp_null(xp);
		fp_null(yp);
		fp18_null(x_q1);
		fp18_null(y_q1);
		fp18_null(x_mq);
		fp18_null(y_mq);
		dv18_null(_t0);
		dv18_null(_t1);
		for(i=0;i<8;i++){
			fp18_null(V1[i]);
		}	
/*========================================================*/
		bn_new(l);
		fp18_new(k0);
		fp18_new(k1);	
		fp18_new(xq);
		fp18_new(yq);
		fp18_new(w_n1);
		fp18_new(w_2n);
		fp18_new(w_21);
		fp18_new(w_11);
		fp18_new(B0);
		for(i=0;i<7;i++)	
			fp18_new(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_new(V2[i]);
		}
		fp3_new(a);
		fp3_new(d);
		fp_new(xp);
		fp_new(yp);
		fp18_new(x_q1);
		fp18_new(y_q1);
		fp18_new(x_mq);
		fp18_new(y_mq);
		dv18_new(_t0);
		dv18_new(_t1);
		for(i=0;i<8;i++){
			fp18_new(V1[i]);
		}
/*========================================================*/
		fp18_set_dig(F,1);
		if (ep_is_infty(P) || ep2_is_infty(Q)) return;
		fp3_copy(a,Q->x);//a=xq
		fp3_copy(d,Q->y);//d=yq
		fp_copy(xp,P->x);//xp
		fp_copy(yp,P->y);//yp

		fp18_zero(xq);
		fp18_zero(yq);
		fp_copy(xq[2][0][0],a[0]);
		fp_copy(xq[2][2][0],a[1]);
		fp_copy(xq[2][1][1],a[2]);
		fp_copy(yq[0][1][0],d[0]);
		fp_copy(yq[0][0][1],d[1]);
		fp_copy(yq[0][2][1],d[2]);	

		fp18_zero(B0);
		fp_copy(B0[0][0][0],B);
		if(bn_sign(m)==1){
			bn_neg(l,m);
			sign=1;	
		}
		else{
			bn_copy(l,m);
		}

/*========================================================*/
/*                      Init data                         */
/*========================================================*/

	/* fill out the first vector of the block */

		fp18_set_dig(V1[3],1);//V1[3]=1
		fp6_dbl(V1[4][0],yq[0]);//the origin w(2,0)=2*yq
		fp6_zero(V1[4][1]);
		fp6_zero(V1[4][2]);
		/* the origin w(3,0)=3*xq^4+12*B*xq=3*xq*(xq^3+4*B) */
		fp6_zero(V1[5][0]);
		fp6_zero(V1[5][1]);		
		fp6_dbl(V1[5][2],xq[2]);
		fp6_add(V1[5][2],V1[5][2],xq[2]);//3*xq
		fp18_sqr(k0,xq);//t0=xq^2
		fp18_mul(k0,k0,xq);//xq^3
		fp6_dbl(k1[0],B0[0]);
		fp6_add(k1[0],k1[0],k1[0]);//k1=4B
		fp6_zero(k1[1]);
		fp6_zero(k1[2]);
		fp18_add(V1[6],k0,k1);
		fp18_mul(V1[5],V1[5],V1[6]);
		/*the origin w(4,0)=4*yq*(xq^6+20Bxq^3-8*B^2))*/
		fp18_sqrn_low(_t0,k0);//xq^6
		fp18_dbl(V1[6],k1);//8B
		fp18_muln_low(_t1,V1[6],B0);//8B^2
		fp18_subc_low(_t0,_t0,_t1);
		fp18_dbl(V1[6],V1[6]);//16B
		fp18_add(k1,k1,V1[6]);
		fp18_muln_low(_t1,k1,k0);
		fp18_addc_low(_t0,_t0,_t1);
		fp18_rdc(V1[6],_t0);
		fp18_dbl(V1[2],V1[4]);//4yq
		fp18_mul(V1[6],V1[6],V1[2]);
		/*the orginal W(5,0)*/
		fp18_sqr(V1[7],V1[5]);
		fp18_mul(V1[7],V1[5],V1[7]);
		fp18_sqr(V1[2],V1[4]);
		fp18_mul(V1[2],V1[2],V1[4]);
		fp18_mul(V1[2],V1[2],V1[6]);
		fp18_sub(V1[7],V1[2],V1[7]);
		fp18_zero(V1[2]);//the orginal W(0,0)=0
		fp18_neg(V1[1],V1[3]);//the orginal W(-1,0)
		fp18_neg(V1[0],V1[4]);//the orginal W(-2,0)

	/* fill out the second vector of the block */

		fp18_set_dig(V2[0],1);// the original W(0,1)
		fp18_set_dig(V2[1],1);// the original W(1,1)
		/* the original W(2,1) =2*xq + xp - ((yp - yq)/(xp-xq))^2*/
		fp18_zero(k0);
		fp6_neg(k0[2],xq[2]);
		fp_copy(k0[0][0][0],xp);
		fp18_inv(k0,k0);//1/(xp-xq)
		fp18_neg(k1,yq);
		fp_copy(k1[0][0][0],yp);
		fp18_mul(k1,k1,k0);
		fp18_sqr(k1,k1);

		fp6_dbl(k0[2],xq[2]);
		fp6_zero(k0[0]);
		fp6_zero(k0[1]);
		fp_copy(k0[0][0][0],xp);//k0=2*xq + xp
		fp18_sub(V2[2],k0,k1);

		fp18_copy(w_n1,xq);
		fp_neg(w_n1[0][0][0],xp);//w_n1=xq-xp
		fp18_sqr(k1,w_n1);//k1=w_n1^2
		fp18_mul(k0,k1,k0);
		fp18_neg(w_2n,k0);
		fp18_copy(k0,yq);
		fp_copy(k0[0][0][0],yp);
		fp18_sqr(k0,k0);//(yp+y_q2)^2
		fp18_add(w_2n,w_2n,k0);	
		fp18_copy(w_21,V2[2]);//the origin w(2,1)
		
		//V2[0]=V2[0]
		fp18_copy(V2[1],w_n1);//V2[1]=V2[1]*w_n1=1*w_n1=w_n1,new w(1,1)
		fp18_mul(V2[2],V2[2],k1);//new w(2,1)

		fp18_inv(k0,k1);//1/w_n1^2
		fp18_mul(w_2n,w_2n,k0);//w(-2,1)
		fp18_set_dig(w_n1,1);//w(1,-1)
		fp18_copy(w_11,V2[1]);//w(1,1)
		
		fp18_inv(initial_data[0],V1[4]);//w(2,0)^(-1)
		fp18_inv(initial_data[1],w_n1);//w_n1=1,w(-1,1)^(-1)
		fp18_inv(initial_data[2],w_2n);//w(-2,1)^(-1)
		fp18_inv(initial_data[3],w_21);//w(2,1)^(-1)
		fp18_inv(initial_data[4],w_11);//w(1,1)^(-1)
		fp18_sqr(initial_data[5],V1[4]);
		fp18_mul(initial_data[6],V1[3],V1[5]);
					
		switch(version){
			case 1:
				precomp(F,x_q1,y_q1,xq,yq,initial_data,V2,V1);
				DoubleAndAdd(initial_data,V2,V1,l);
				break;
			case 2:
				precomp_lz(F,x_q1,y_q1,xq,yq,initial_data,V2,V1);
				DoubleAndAdd_lz(initial_data,V2,V1,l);
				break;
			default:
				precomp(F,x_q1,y_q1,xq,yq,initial_data,V2,V1);
				DoubleAndAdd(initial_data,V2,V1,l);
				break;
		}
	/*x_mq = xq - (V1[2] * V1[4]) / (V1[3] * V1[3])*/
		fp18_sqr(k0,V1[3]);
		fp18_inv(x_mq,k0);
		fp18_mul(k1,V1[2],V1[4]);	
		fp18_mul(k1,k1,x_mq);
		fp18_sub(x_mq,xq,k1);//notice that when m=r,V1[3]=0
		
	/*y_mq=(V1[2]^2*V1[5]-V1[4]^2*V1[1])/4*yq*V1[3]^3*/
		fp18_mul(k0,k0,V1[3]);
		fp6_zero(k1[1]);
		fp6_zero(k1[2]);
		fp6_dbl(k1[0],yq[0]);
		fp6_dbl(k1[0],k1[0]);
		fp18_mul(y_mq,k1,k0);
		fp18_inv(y_mq,y_mq);
		fp18_sqr(k0,V1[2]);
		fp18_mul(k0,k0,V1[5]);
		fp18_sqr(k1,V1[4]);
		fp18_mul(k1,k1,V1[1]);
		fp18_sub(k0,k0,k1);
		fp18_mul(y_mq,k0,y_mq);

		/* f_{-a,Q}(P) = 1/f_{a,Q}(P). */
		if(sign==1){
			fp18_inv(k0, V2[1]);
			fp18_neg(y_mq,y_mq);
		}
		else{
			fp18_copy(k0,V2[1]);	
		}

		fp18_frb(F,F,1);//l^p
		fp18_frb(x_q1,x_q1,1);//x_q1^p
		fp18_frb(y_q1,y_q1,1);//y_q1^p
		fp18_mul(F,F,k0);
		
		fp18_zero(k0);
		fp_sub(k0[2][0][0],x_mq[2][0][0], x_q1[2][0][0]);
		fp_sub(k0[2][2][0],x_mq[2][2][0], x_q1[2][2][0]);
		fp_sub(k0[2][1][1],x_mq[2][1][1], x_q1[2][1][1]);	
		fp18_inv(k0,k0);
		fp18_copy(k1,y_q1);
		fp_sub(k1[0][1][0],y_mq[0][1][0], y_q1[0][1][0]);
		fp_sub(k1[0][0][1],y_mq[0][0][1], y_q1[0][0][1]);
		fp_sub(k1[0][2][1],y_mq[0][2][1], y_q1[0][2][1]);	
		fp18_mul_dxs(k0,k0,k1);
		fp_sub(x_q1[0][0][0],x_q1[0][0][0],xp);
		fp18_mul(k0,k0,x_q1);
		fp18_neg(k0,k0);
		fp_sub(y_q1[0][0][0],y_q1[0][0][0],yp);
		fp18_add(k0,k0,y_q1);
		fp18_mul(F,F,k0);
		
		pp_exp_k18(F,F);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}	

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {
		bn_free(l);
		fp18_free(k0);
		fp18_free(k1);	
		fp18_free(xq);
		fp18_free(yq);
		fp18_free(w_n1);
		fp18_free(w_2n);
		fp18_free(w_21);
		fp18_free(w_11);
		fp18_free(B0);
		for(i=0;i<7;i++)	
			fp18_free(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_free(V2[i]);
		}
		fp3_free(a);
		fp3_free(d);
		fp_free(xp);
		fp_free(yp);
		fp18_free(x_q1);
		fp18_free(y_q1);
		fp18_free(x_mq);
		fp18_free(y_mq);
		dv18_free(_t0);
		dv18_free(_t1);
		for(i=0;i<8;i++){
			fp18_free(V1[i]);
		}
	}
}
