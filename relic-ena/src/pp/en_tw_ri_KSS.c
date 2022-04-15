/*==========================================================================*/
/*==========                                                   =============*/
/*          Improved Elliptic Net Alg on Twisted KSS curves (RI)            */
/*==========                                                   =============*/
/*==========================================================================*/


#include"ena.h"

static void fp18_mul_fp3(fp18_t c,fp3_t b,fp18_t a){
	fp6_t t0;
	fp6_null(t0);
	fp6_new(t0);
	
	fp6_zero(t0);
	fp_copy(t0[0][0],b[0]);
	fp_copy(t0[1][1],b[2]);
	fp_copy(t0[2][0],b[1]);

	fp6_mul(c[0],a[0],t0);
	fp6_mul(c[1],a[1],t0);
	fp6_mul(c[2],a[2],t0);

	fp6_free(t0);
}
static void fp18_muln_fp3(dv18_t c,fp3_t b,fp18_t a){
	fp6_t t0;
	fp6_null(t0);
	fp6_new(t0);
	
	fp6_zero(t0);
	fp_copy(t0[0][0],b[0]);
	fp_copy(t0[1][1],b[2]);
	fp_copy(t0[2][0],b[1]);

	fp6_muln_low(c[0],a[0],t0);
	fp6_muln_low(c[1],a[1],t0);
	fp6_muln_low(c[2],a[2],t0);

	fp6_free(t0);
}

static void precomp(fp18_t F,fp18_t x_q1,fp18_t y_q1,fp3_t xq,fp3_t yq,fp18_t initial_data[],fp18_t V2[],fp3_t V1[]){
	fp18_t beta,gamma1,S_0,P_0,t0,tv2[3];
    	fp3_t alpha,_S1[5],_P1[5],tv1[7],k0,w2,w13,vt1,vt2,vt3;
	int i;

/*========================================================*/
/*                          Init                          */
/*========================================================*/
	TRY{
		fp3_null(alpha);
		fp3_null(k0);
		fp3_null(w2);
		fp3_null(w13);
		fp3_null(vt1);
		fp3_null(vt2);
		fp3_null(vt3);
		fp18_null(beta);
		fp18_null(gamma1);
		fp18_null(S_0);
		fp18_null(P_0);
		fp18_null(t0);

		for(i=0;i<5;i++){
			fp3_null(_S1[i]);
			fp3_null(_P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_null(tv2[i]);
		}
		for(i=0;i<7;i++){
			fp3_null(tv1[i]);
		}
/*========================================================*/
		fp3_new(alpha);
		fp3_new(k0);
		fp18_new(beta);
		fp18_new(gamma1);
		fp18_new(S_0);
		fp18_new(P_0);
		fp18_new(t0);
		fp3_new(w2);
		fp3_new(w13);
		fp3_new(vt1);
		fp3_new(vt2);
		fp3_new(vt3);

		for(i=0;i<5;i++){
			fp3_new(_S1[i]);
			fp3_new(_P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_new(tv2[i]);
		}
		for(i=0;i<7;i++){
			fp3_new(tv1[i]);
		}
/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

		for(i=0;i<3;i++){
			fp18_copy(tv2[i],V2[i]);
		}
		for(i=0;i<7;i++){
			fp3_copy(tv1[i],V1[i]);
		}
	
		//initial data contains the precomputed inverses
		fp_copy(alpha[0],initial_data[0][0][0][0]);
		fp_copy(alpha[1],initial_data[0][0][2][0]);
		fp_copy(alpha[2],initial_data[0][0][1][1]);
		fp18_copy(beta,initial_data[1]);// inverse of W(-1,1)
		fp18_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
		fp_copy(w2[0],initial_data[5][0][0][0]);// inverse of W(1,1)
		fp_copy(w2[1],initial_data[5][0][2][0]);
		fp_copy(w2[2],initial_data[5][0][1][1]);
		fp_copy(w13[0],initial_data[6][0][0][0]);// inverse of W(1,1)
		fp_copy(w13[1],initial_data[6][0][2][0]);
		fp_copy(w13[2],initial_data[6][0][1][1]);

		fp18_sqr(S_0,tv2[1]);//1S_2
		fp18_mul(P_0,tv2[0],tv2[2]);//1M_2
	    
		for (int j = 0; j <= 4; j++)//6S,6M
		{
			fp3_sqr(_S1[j],tv1[j + 1]);
			fp3_mul(_P1[j],tv1[j],tv1[j + 2]);
		}
			
		for (int j = 1; j <= 3; j++)
		{
		/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
			fp3_mul(tv1[2*j-2],_S1[j-1],_P1[j+1]);
			fp3_mul(k0,_S1[j+1],_P1[j-1]);
			fp3_sub(tv1[2*j-2],tv1[2*j-2],k0);
			fp3_mul(tv1[2*j-2],tv1[2*j-2],alpha);
		/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
			fp3_mul(tv1[2*j-1],_S1[j],_P1[j+1]);		
			fp3_mul(k0,_S1[j+1],_P1[j]);
			fp3_sub(tv1[2*j-1],tv1[2*j-1],k0);
		}

		/* V1[6] = (vt1*w2 - vt2 * w13) / V1[2] */
		fp3_mul(vt1,tv1[3],tv1[5]);
		fp3_sqr(vt2,tv1[4]);
		fp3_copy(vt3,tv1[2]);
		fp3_mul(tv1[6],vt1,w2);
		fp3_mul(k0,vt2,w13);
		fp3_sub(tv1[6],tv1[6],k0);
		for(int j=0;j<=5;j++)fp3_mul(tv1[j],tv1[j],vt3);

		/* tv2[0] = S1[2] * P_0 - P1[2] * S_0 */
		fp18_mul_fp3(tv2[0],_S1[2],P_0);
		fp18_mul_fp3(t0,_P1[2],S_0);
		fp18_sub(tv2[0],tv2[0],t0);
		/* tv2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
		fp18_mul_fp3(tv2[1],_S1[3],P_0);
		fp18_mul_fp3(t0,_P1[3],S_0);
		fp18_sub(tv2[1],tv2[1],t0);
		fp18_mul_dxs(tv2[1],tv2[1],beta);
		/* tv2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
		fp18_mul_fp3(tv2[2],_P1[4],S_0);
		fp18_mul_fp3(t0,_S1[4],P_0);
		fp18_sub(tv2[2],tv2[2],t0);
		fp18_mul(tv2[2],tv2[2],gamma1);

	/*x_q1 = xq - (tv1[2] * tv1[4]) / (tv1[3] * tv1[3])*/
		fp3_sqr(alpha,tv1[3]);
		fp3_inv(_S1[0],alpha);
		fp3_mul(_S1[1],tv1[2],tv1[4]);	
		fp3_mul(_S1[0],_S1[0],_S1[1]);
		fp3_sub(_S1[0],xq,_S1[0]);
		//notice that when m=r,tv1[3]=0	
		fp18_zero(x_q1);
		fp_copy(x_q1[2][0][0],_S1[0][0]);
		fp_copy(x_q1[2][2][0],_S1[0][1]);
		fp_copy(x_q1[2][1][1],_S1[0][2]);

	/*y_q1=(tv1[2]^2*tv1[5]-tv1[4]^2*tv1[1])/4*yq*tv1[3]^3*/
		fp3_mul(alpha,alpha,tv1[3]);
		fp3_dbl(_S1[0],yq);
		fp3_dbl(_S1[0],_S1[0]);
		fp3_mul(_S1[0],_S1[0],alpha);
		fp3_inv(alpha,_S1[0]);
		fp3_sqr(_S1[0],tv1[2]);
		fp3_mul(_S1[0],_S1[0],tv1[5]);
		fp3_sqr(_S1[1],tv1[4]);
		fp3_mul(_S1[1],_S1[1],tv1[1]);
		fp3_sub(_S1[0],_S1[0],_S1[1]);
		fp3_mul(_S1[0],_S1[0],alpha);
		fp18_zero(y_q1);
		fp_copy(y_q1[0][1][0],_S1[0][0]);
		fp_copy(y_q1[0][0][1],_S1[0][1]);
		fp_copy(y_q1[0][2][1],_S1[0][2]);

		fp18_copy(F,tv2[1]);
}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {
		fp3_free(alpha);	
		fp3_free(k0);	
		fp18_free(beta);
		fp18_free(gamma1);
		fp18_free(S_0);
		fp18_free(P_0);
		fp18_free(t0);
		fp3_free(w2);
		fp3_free(w13);
		fp3_free(vt1);
		fp3_free(vt2);
		fp3_free(vt3);

		for(i=0;i<6;i++){
			fp3_free(_S1[i]);
			fp3_free(_P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_free(tv2[i]);
		}
		for(i=0;i<8;i++){
			fp3_free(tv1[i]);
		}
	}
}


static void precomp_lz(fp18_t F,fp18_t x_q1,fp18_t y_q1,fp3_t xq,fp3_t yq,fp18_t initial_data[],fp18_t V2[],fp3_t V1[]){
	fp18_t beta,gamma1,S_0,P_0,tv2[3];
    	fp3_t alpha,_S1[5],_P1[5],tv1[7],w2,w13,vt1,vt2,vt3;
	int i;
	dv18_t t0,t1;

/*========================================================*/
/*                          Init                          */
/*========================================================*/
	TRY{
		fp3_null(alpha);
		fp18_null(beta);
		fp18_null(gamma1);
		fp18_null(S_0);
		fp18_null(P_0);
		dv18_null(t0);
		dv18_null(t1);
		fp3_null(w2);
		fp3_null(w13);
		fp3_null(vt1);
		fp3_null(vt2);
		fp3_null(vt3);

		for(i=0;i<5;i++){
			fp3_null(_S1[i]);
			fp3_null(_P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_null(tv2[i]);
		}
		for(i=0;i<7;i++){
			fp3_null(tv1[i]);
		}
/*========================================================*/
		fp3_new(alpha);
		fp18_new(beta);
		fp18_new(gamma1);
		fp18_new(S_0);
		fp18_new(P_0);
		dv18_new(t0);
		dv18_new(t1);
		fp3_new(w2);
		fp3_new(w13);
		fp3_new(vt1);
		fp3_new(vt2);
		fp3_new(vt3);


		for(i=0;i<5;i++){
			fp3_new(_S1[i]);
			fp3_new(_P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_new(tv2[i]);
		}
		for(i=0;i<7;i++){
			fp3_new(tv1[i]);
		}
/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

		for(i=0;i<3;i++){
			fp18_copy(tv2[i],V2[i]);
		}
		for(i=0;i<7;i++){
			fp3_copy(tv1[i],V1[i]);
		}
	
		//initial data contains the precomputed inverses
		fp_copy(alpha[0],initial_data[0][0][0][0]);
		fp_copy(alpha[1],initial_data[0][0][2][0]);
		fp_copy(alpha[2],initial_data[0][0][1][1]);
		fp18_copy(beta,initial_data[1]);// inverse of W(-1,1)
		fp18_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
		fp_copy(w2[0],initial_data[5][0][0][0]);// inverse of W(1,1)
		fp_copy(w2[1],initial_data[5][0][2][0]);
		fp_copy(w2[2],initial_data[5][0][1][1]);
		fp_copy(w13[0],initial_data[6][0][0][0]);// inverse of W(1,1)
		fp_copy(w13[1],initial_data[6][0][2][0]);
		fp_copy(w13[2],initial_data[6][0][1][1]);

		fp18_sqr(S_0,tv2[1]);//1S_2
		fp18_mul(P_0,tv2[0],tv2[2]);//1M_2
	    
		for (int j = 0; j <= 4; j++)//6S,6M
		{
			fp3_sqr(_S1[j],tv1[j + 1]);
			fp3_mul(_P1[j],tv1[j],tv1[j + 2]);
		}
			
		for (int j = 1; j <= 3; j++)
		{
		/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
			fp3_muln_low(t0[0][0],_S1[j-1],_P1[j+1]);
			fp3_muln_low(t1[0][0],_S1[j+1],_P1[j-1]);
			fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
			fp3_rdc(tv1[2*j-2],t0[0][0]);
			fp3_mul(tv1[2*j-2],tv1[2*j-2],alpha);
		/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
			fp3_muln_low(t0[0][0],_S1[j],_P1[j+1]);		
			fp3_muln_low(t1[0][0],_S1[j+1],_P1[j]);
			fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
			fp3_rdc(tv1[2*j-1],t0[0][0]);
		}

		/* V1[6] = (vt1*w2 - vt2 * w13) / V1[2] */

		fp3_mul(vt1,tv1[3],tv1[5]);
		fp3_sqr(vt2,tv1[4]);
		fp3_copy(vt3,tv1[2]);
		fp3_muln_low(t0[0][0],vt1,w2);
		fp3_muln_low(t1[0][0],vt2,w13);
		fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
		fp3_rdc(tv1[6],t0[0][0]);
		for(int j=0;j<=5;j++)fp3_mul(tv1[j],tv1[j],vt3);

		/* tv2[0] = S1[2] * P_0 - P1[2] * S_0 */
		fp18_muln_fp3(t0,_S1[2],P_0);
		fp18_muln_fp3(t1,_P1[2],S_0);
		fp18_subc_low(t0,t0,t1);
		fp18_rdc(tv2[0],t0);
		/* tv2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
		fp18_muln_fp3(t0,_S1[3],P_0);
		fp18_muln_fp3(t1,_P1[3],S_0);
		fp18_subc_low(t0,t0,t1);
		fp18_rdc(tv2[1],t0);
		fp18_mul_dxs(tv2[1],tv2[1],beta);
		/* tv2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
		fp18_muln_fp3(t0,_P1[4],S_0);
		fp18_muln_fp3(t1,_S1[4],P_0);
		fp18_subc_low(t0,t0,t1);
		fp18_rdc(tv2[2],t0);
		fp18_mul(tv2[2],tv2[2],gamma1);

	/*x_q1 = xq - (tv1[2] * tv1[4]) / (tv1[3] * tv1[3])*/
		fp3_sqr(alpha,tv1[3]);
		fp3_inv(_S1[0],alpha);
		fp3_mul(_S1[1],tv1[2],tv1[4]);	
		fp3_mul(_S1[0],_S1[0],_S1[1]);
		fp3_sub(_S1[0],xq,_S1[0]);
		//notice that when m=r,tv1[3]=0	
		fp18_zero(x_q1);
		fp_copy(x_q1[2][0][0],_S1[0][0]);
		fp_copy(x_q1[2][2][0],_S1[0][1]);
		fp_copy(x_q1[2][1][1],_S1[0][2]);

	/*y_q1=(tv1[2]^2*tv1[5]-tv1[4]^2*tv1[1])/4*yq*tv1[3]^3*/
		fp3_mul(alpha,alpha,tv1[3]);
		fp3_dbl(_S1[0],yq);
		fp3_dbl(_S1[0],_S1[0]);
		fp3_mul(_S1[0],_S1[0],alpha);
		fp3_inv(alpha,_S1[0]);
		fp3_sqr(_S1[0],tv1[2]);
		fp3_mul(_S1[0],_S1[0],tv1[5]);
		fp3_sqr(_S1[1],tv1[4]);
		fp3_mul(_S1[1],_S1[1],tv1[1]);
		fp3_sub(_S1[0],_S1[0],_S1[1]);
		fp3_mul(_S1[0],_S1[0],alpha);
		fp18_zero(y_q1);
		fp_copy(y_q1[0][1][0],_S1[0][0]);
		fp_copy(y_q1[0][0][1],_S1[0][1]);
		fp_copy(y_q1[0][2][1],_S1[0][2]);

		fp18_copy(F,tv2[1]);
}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {
		fp3_free(alpha);	
		fp18_free(beta);
		fp18_free(gamma1);
		fp18_free(S_0);
		fp18_free(P_0);
		dv18_free(t0);
		dv18_free(t1);
		fp3_free(w2);
		fp3_free(w13);
		fp3_free(vt1);
		fp3_free(vt2);
		fp3_free(vt3);

		for(i=0;i<6;i++){
			fp3_free(_S1[i]);
			fp3_free(_P1[i]);
		}
		for(i=0;i<3;i++){
			fp18_free(tv2[i]);
		}
		for(i=0;i<8;i++){
			fp3_free(tv1[i]);
		}
	}
}

/*double and add alg with 7 vector*/
static void DoubleAndAdd_2(fp18_t initial_data[],fp18_t V2[],fp3_t V1[],bn_t m){
	fp18_t beta,gamma1,gamma2,delta,S_0,P_0,t0;
    	fp3_t alpha,w2,w13,vt1,vt2,vt3,_S1[5],_P1[5];
	int nb,i,len = bn_bits(m) + 1;
	int8_t* naf=RLC_ALLOCA(int8_t, len);
	bn_rec_naf(naf,&len, m, 2);

/*========================================================*/
/*                          Init                          */
/*========================================================*/
	TRY {
		fp3_null(alpha);
		fp3_null(vt1);
		fp3_null(vt2);
		fp3_null(vt3);
		fp3_null(w2);
		fp3_null(w13);
		fp18_null(beta);
		fp18_null(gamma1);
		fp18_null(gamma2);
		fp18_null(delta);
		fp18_null(S_0);
		fp18_null(P_0);
		fp18_null(t0);
		for(i=0;i<5;i++){
			fp3_null(_S1[i]);
			fp3_null(_P1[i]);
		}
/*========================================================*/
		fp3_new(alpha);
		fp3_new(vt1);
		fp3_new(vt2);
		fp3_new(vt3);
		fp3_new(w2);
		fp3_new(w13);
		fp18_new(beta);
		fp18_new(gamma1);
		fp18_new(gamma2);
		fp18_new(delta);
		fp18_new(S_0);
		fp18_new(P_0);
		fp18_new(t0);
	
		for(i=0;i<5;i++){
			fp3_new(_S1[i]);
			fp3_new(_P1[i]);
		}

/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

		nb = len - 1;
		//initial data contains the precomputed inverses
		fp_copy(alpha[0],initial_data[0][0][0][0]);
		fp_copy(alpha[1],initial_data[0][0][2][0]);
		fp_copy(alpha[2],initial_data[0][0][1][1]);
		fp18_copy(beta,initial_data[1]);// inverse of W(-1,1)
		fp18_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
		fp18_copy(gamma2,initial_data[3]);// inverse of W(-2,1)
		fp18_copy(delta,initial_data[4]); // inverse of W(1,1)
		fp_copy(w2[0],initial_data[5][0][0][0]);// inverse of W(1,1)
		fp_copy(w2[1],initial_data[5][0][2][0]);
		fp_copy(w2[2],initial_data[5][0][1][1]);
		fp_copy(w13[0],initial_data[6][0][0][0]);// inverse of W(1,1)
		fp_copy(w13[1],initial_data[6][0][2][0]);
		fp_copy(w13[2],initial_data[6][0][1][1]);

		for(i=nb-1;i>=0;i--){
			int b=naf[i];	
			fp18_sqr(S_0,V2[1]);//1S_2
			fp18_mul(P_0,V2[0],V2[2]);//1M_2
			for (int j = 0; j <= 4; j++)//5S,5M
			{
                		fp3_sqr(_S1[j],V1[j + 1]);
                		fp3_mul(_P1[j],V1[j],V1[j + 2]);
			}
			
			if (b == 0)
			{
                	//update i to 2i
				for (int j = 1; j <= 3; j++)//j=3 or 4
				{
				/* S[j - 1] * P1[j] - S[j] * P1[j - 1] */
					fp3_mul(V1[2*j-2],_S1[j-1],_P1[j]);
					fp3_mul(t0[0][0],_S1[j],_P1[j-1]);
					fp3_sub(V1[2*j-2],V1[2*j-2],t0[0][0]);
				/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
					fp3_mul(V1[2*j-1],_S1[j-1],_P1[j+1]);	
					fp3_mul(t0[0][0],_S1[j+1],_P1[j-1]);
					fp3_sub(V1[2*j-1],V1[2*j-1],t0[0][0]);
					fp3_mul(V1[2*j-1],V1[2*j-1],alpha);
				}
	
				/* V1[6] = S[3] * P1[4] - S[4] * P1[3] */
				fp3_mul(V1[6],_S1[3],_P1[4]);
				fp3_mul(t0[0][0],_S1[4],_P1[3]);
				fp3_sub(V1[6],V1[6],t0[0][0]);

				/* V2[0] = (S1[1] * P_0 - P1[1] * S_0) * delta*/
				fp18_mul_fp3(V2[0],_S1[1],P_0);
				fp18_mul_fp3(t0,_P1[1],S_0);

				fp18_sub(V2[0],V2[0],t0);
				fp18_mul(V2[0],V2[0],delta);
				/*V2[1] = S1[2] * P_0 - P1[2] * S_0*/
				fp18_mul_fp3(V2[1],_S1[2],P_0);
				fp18_mul_fp3(t0,_P1[2],S_0);
				fp18_sub(V2[1],V2[1],t0);
				/* V2[2] = (S1[3] * P_0 - P1[3] * S_0) * beta */
				fp18_mul_fp3(V2[2],_S1[3],P_0);
				fp18_mul_fp3(t0,_P1[3],S_0);
				fp18_sub(V2[2],V2[2],t0);
				fp18_mul_dxs(V2[2],V2[2],beta);
			}
			else if(b==-1)
{
				for (int j = 1; j <= 3; j++)
				{
				/* (S[j - 1] * P1[j] - S[j] * P1[j - 1]) */
					fp3_mul(V1[2*j-1],_S1[j-1],_P1[j]);
					fp3_mul(t0[0][0],_S1[j],_P1[j-1]);
					fp3_sub(V1[2*j-1],V1[2*j-1],t0[0][0]);
				/* (S[j-1] * P1[j + 1] - S[j + 1] * P1[j-1])*alpha */
					fp3_mul(V1[2*j],_S1[j-1],_P1[j+1]);		
					fp3_mul(t0[0][0],_S1[j+1],_P1[j-1]);
					fp3_sub(V1[2*j],V1[2*j],t0[0][0]);
					fp3_mul(V1[2*j],V1[2*j],alpha);

				}
				/* V1[0] = (vt1*w2 - vt2 * w13) / V1[4] */
				fp3_mul(vt1,V1[3],V1[1]);
				fp3_sqr(vt2,V1[2]);
				fp3_copy(vt3,V1[4]);
				fp3_mul(V1[0],vt1,w2);
				fp3_mul(t0[0][0],vt2,w13);
				fp3_sub(V1[0],V1[0],t0[0][0]);
				for(int j=1;j<=6;j++)fp3_mul(V1[j],V1[j],vt3);

				/* V2[0] = (S1[0] * P_0 - P1[0] * S_0 )*gamma2*/
				fp18_mul_fp3(V2[0],_S1[0],P_0);
				fp18_mul_fp3(t0,_P1[0],S_0);
				fp18_sub(V2[0],V2[0],t0);
				fp18_mul(V2[0],V2[0],gamma2);
				/* V2[1] = (S1[1] * P_0 - P1[1] * S_0) * delta */
				fp18_mul_fp3(V2[1],_S1[1],P_0);
				fp18_mul_fp3(t0,_P1[1],S_0);
				fp18_sub(V2[1],V2[1],t0);
				fp18_mul(V2[1],V2[1],delta);
				/* V2[2] = (S1[2] * P_0 - P1[2] * S_0) */
				fp18_mul_fp3(V2[2],_S1[2],P_0);
				fp18_mul_fp3(t0,_P1[2],S_0);
				fp18_sub(V2[2],V2[2],t0);
			}
			else
			{
				for (int j = 1; j <= 3; j++)
				{
				/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
					fp3_mul(V1[2*j-2],_S1[j-1],_P1[j+1]);
					fp3_mul(t0[0][0],_S1[j+1],_P1[j-1]);
					fp3_sub(V1[2*j-2],V1[2*j-2],t0[0][0]);
					fp3_mul(V1[2*j-2],V1[2*j-2],alpha);
				/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
					fp3_mul(V1[2*j-1],_S1[j],_P1[j+1]);		
					fp3_mul(t0[0][0],_S1[j+1],_P1[j]);
					fp3_sub(V1[2*j-1],V1[2*j-1],t0[0][0]);
				}
				//cout << "V1[3]=" << V1[3] << endl;
				//cout << "V1[5]=" << V1[5] << endl;
				/* V1[6] = (vt1*w2 - vt2 * w13) / V1[2] */
				fp3_mul(vt1,V1[3],V1[5]);
				fp3_sqr(vt2,V1[4]);
				fp3_copy(vt3,V1[2]);
				fp3_mul(V1[6],vt1,w2);
				fp3_mul(t0[0][0],vt2,w13);
				fp3_sub(V1[6],V1[6],t0[0][0]);
				for(int j=0;j<=5;j++)fp3_mul(V1[j],V1[j],vt3);

				/* V2[0] = S1[2] * P_0 - P1[2] * S_0 */
				fp18_mul_fp3(V2[0],_S1[2],P_0);
				fp18_mul_fp3(t0,_P1[2],S_0);
				fp18_sub(V2[0],V2[0],t0);
				/* V2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
				fp18_mul_fp3(V2[1],_S1[3],P_0);
				fp18_mul_fp3(t0,_P1[3],S_0);
				fp18_sub(V2[1],V2[1],t0);
				fp18_mul_dxs(V2[1],V2[1],beta);
				/* V2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
				fp18_mul_fp3(V2[2],_P1[4],S_0);
				fp18_mul_fp3(t0,_S1[4],P_0);
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
		fp3_free(alpha);	
		fp3_free(vt1);
		fp3_free(vt2);
		fp3_free(vt3);
		fp3_free(w2);
		fp3_free(w13);
		fp18_free(beta);
		fp18_free(gamma1);
		fp18_free(gamma2);
		fp18_free(delta);
		fp18_free(S_0);
		fp18_free(P_0);
		fp18_free(t0);
	
		for(i=0;i<5;i++){
			fp3_free(_S1[i]);
			fp3_free(_P1[i]);
		}
}
}

/*improve double and add alg with lazy reduction*/
static void DoubleAndAdd_2_lz(fp18_t initial_data[],fp18_t V2[],fp3_t V1[],bn_t m){
    	fp18_t beta,gamma1,gamma2,delta,S_0,P_0;
	fp3_t alpha,w2,w13,vt1,vt2,vt3,_S1[5],_P1[5];
	dv18_t t0,t1;
	int nb,i,len = bn_bits(m) + 1;
	int8_t* naf=RLC_ALLOCA(int8_t, len);
	bn_rec_naf(naf,&len, m, 2);

/*========================================================*/
/*                        Init                            */
/*========================================================*/
	TRY{
		fp3_null(alpha);
		fp3_null(vt1);
		fp3_null(vt2);
		fp3_null(vt3);
		fp3_null(w2);
		fp3_null(w13);
		fp18_null(beta);
		fp18_null(gamma1);
		fp18_null(gamma2);
		fp18_null(delta);
		fp18_null(S_0);
		fp18_null(P_0);
		dv18_null(t0);
		dv18_null(t1);
		for(i=0;i<5;i++){
			fp3_null(_S1[i]);
			fp3_null(_P1[i]);
		}
/*========================================================*/
		for(i=0;i<5;i++){
			fp3_null(_S1[i]);
			fp3_null(_P1[i]);
		}
		fp3_new(alpha);
		fp3_new(vt1);
		fp3_new(vt2);
		fp3_new(vt3);
		fp3_new(w2);
		fp3_new(w13);
		fp18_new(beta);
		fp18_new(gamma1);
		fp18_new(gamma2);
		fp18_new(delta);
		fp18_new(S_0);
		fp18_new(P_0);
	
		for(i=0;i<5;i++){
			fp3_new(_S1[i]);
			fp3_new(_P1[i]);
		}
        
		dv18_new(t0);
		dv18_new(t1);
	
/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

		nb = len - 1;
		//initial data contains the precomputed inverses
		fp_copy(alpha[0],initial_data[0][0][0][0]);
		fp_copy(alpha[1],initial_data[0][0][2][0]);
		fp_copy(alpha[2],initial_data[0][0][1][1]);
		fp18_copy(beta,initial_data[1]);// inverse of W(-1,1)
		fp18_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
		fp18_copy(gamma2,initial_data[3]);// inverse of W(-2,1)
		fp18_copy(delta,initial_data[4]); // inverse of W(1,1)
		fp_copy(w2[0],initial_data[5][0][0][0]);// inverse of W(1,1)
		fp_copy(w2[1],initial_data[5][0][2][0]);
		fp_copy(w2[2],initial_data[5][0][1][1]);
		fp_copy(w13[0],initial_data[6][0][0][0]);// inverse of W(1,1)
		fp_copy(w13[1],initial_data[6][0][2][0]);
		fp_copy(w13[2],initial_data[6][0][1][1]);
		
		for(i=nb-1;i>=0;i--){
			int b=naf[i];	
			fp18_sqr(S_0,V2[1]);//1S_2
			fp18_mul(P_0,V2[0],V2[2]);//1M_2
			
			for (int j = 0; j <= 4; j++)//5S,5M
			{
                		fp3_sqr(_S1[j],V1[j + 1]);
                		fp3_mul(_P1[j],V1[j],V1[j + 2]);
			}
			
			if (b == 0)
			{
				for (int j = 1; j <= 3; j++)//j=3 or 4
				{
				/* S[j - 1] * P1[j] - S[j] * P1[j - 1] */
					fp3_muln_low(t0[0][0],_S1[j-1],_P1[j]);
					fp3_muln_low(t1[0][0],_S1[j],_P1[j-1]);
					fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
					fp3_rdc(V1[2*j-2],t0[0][0]);
				/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
					fp3_muln_low(t0[0][0],_S1[j-1],_P1[j+1]);	
					fp3_muln_low(t1[0][0],_S1[j+1],_P1[j-1]);
					fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
					fp3_rdc(V1[2*j-1],t0[0][0]);
					fp3_mul(V1[2*j-1],V1[2*j-1],alpha);
				}
	
				/* V1[6] = S[3] * P1[4] - S[4] * P1[3] */
				fp3_muln_low(t0[0][0],_S1[3],_P1[4]);
				fp3_muln_low(t1[0][0],_S1[4],_P1[3]);
				fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
				fp3_rdc(V1[6],t0[0][0]);
	
				/* V2[0] = (S1[1] * P_0 - P1[1] * S_0) * delta*/
				fp18_muln_fp3(t0,_S1[1],P_0);
				fp18_muln_fp3(t1,_P1[1],S_0);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V2[0],t0);
				fp18_mul(V2[0],V2[0],delta);
				/*V2[1] = S1[2] * P_0 - P1[2] * S_0*/
				fp18_muln_fp3(t0,_S1[2],P_0);
				fp18_muln_fp3(t1,_P1[2],S_0);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V2[1],t0);
				/* V2[2] = (S1[3] * P_0 - P1[3] * S_0) * beta */
				fp18_muln_fp3(t0,_S1[3],P_0);
				fp18_muln_fp3(t1,_P1[3],S_0);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V2[2],t0);
				fp18_mul_dxs(V2[2],V2[2],beta);
			}
			else if(b==-1)
			{
				for (int j = 1; j <= 3; j++)
				{
				/* (S[j - 1] * P1[j] - S[j] * P1[j - 1]) */
					fp3_muln_low(t0[0][0],_S1[j-1],_P1[j]);
					fp3_muln_low(t1[0][0],_S1[j],_P1[j-1]);
					fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
					fp3_rdc(V1[2*j-1],t0[0][0]);
				/* S[j-1] * P1[j + 1] - S[j + 1] * P1[j-1] */
					fp3_muln_low(t0[0][0],_S1[j-1],_P1[j+1]);		
					fp3_muln_low(t1[0][0],_S1[j+1],_P1[j-1]);
					fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
					fp3_rdc(V1[2*j],t0[0][0]);
					fp3_mul(V1[2*j],V1[2*j],alpha);
				}
				//cout << "V1[3]=" << V1[3] << endl;
				//cout << "V1[5]=" << V1[5] << endl;
	
				/* V1[0] = (vt1*w2 - vt2 * w13) / V1[4] */

				fp3_mul(vt1,V1[3],V1[1]);
				fp3_sqr(vt2,V1[2]);
				fp3_copy(vt3,V1[4]);
				fp3_muln_low(t0[0][0],vt1,w2);
				fp3_muln_low(t1[0][0],vt2,w13);
				fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
				fp3_rdc(V1[0],t0[0][0]);
				for(int j=1;j<=6;j++)fp3_mul(V1[j],V1[j],vt3);


				/* V2[0] = (S1[0] * P_0 - P1[0] * S_0)* gamma2  */
				fp18_muln_fp3(t0,_S1[0],P_0);
				fp18_muln_fp3(t1,_P1[0],S_0);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V2[0],t0);
				fp18_mul(V2[0],V2[0],gamma2);
				/* V2[1] = (S1[1] * P_0 - P1[1] * S_0) * delta */
				fp18_muln_fp3(t0,_S1[1],P_0);
				fp18_muln_fp3(t1,_P1[1],S_0);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V2[1],t0);
				fp18_mul(V2[1],V2[1],delta);
				/* V2[2] = (S1[2] * P_0 - P1[2] * S_0) */
				fp18_muln_fp3(t0,_S1[2],P_0);
				fp18_muln_fp3(t1,_P1[2],S_0);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V2[2],t0);

			}
			else
			{
				for (int j = 1; j <= 3; j++)
				{
				/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
					fp3_muln_low(t0[0][0],_S1[j-1],_P1[j+1]);
					fp3_muln_low(t1[0][0],_S1[j+1],_P1[j-1]);
					fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
					fp3_rdc(V1[2*j-2],t0[0][0]);
					fp3_mul(V1[2*j-2],V1[2*j-2],alpha);
				/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
					fp3_muln_low(t0[0][0],_S1[j],_P1[j+1]);		
					fp3_muln_low(t1[0][0],_S1[j+1],_P1[j]);
					fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
					fp3_rdc(V1[2*j-1],t0[0][0]);
				}
				//cout << "V1[3]=" << V1[3] << endl;
				//cout << "V1[5]=" << V1[5] << endl;
	
				/* V1[6] = (vt1*w2 - vt2 * w13) / V1[2] */

				fp3_mul(vt1,V1[3],V1[5]);
				fp3_sqr(vt2,V1[4]);
				fp3_copy(vt3,V1[2]);
				fp3_muln_low(t0[0][0],vt1,w2);
				fp3_muln_low(t1[0][0],vt2,w13);
				fp3_subc_low(t0[0][0],t0[0][0],t1[0][0]);
				fp3_rdc(V1[6],t0[0][0]);
				for(int j=0;j<=5;j++)fp3_mul(V1[j],V1[j],vt3);


				/* V2[0] = S1[2] * P_0 - P1[2] * S_0 */
				fp18_muln_fp3(t0,_S1[2],P_0);
				fp18_muln_fp3(t1,_P1[2],S_0);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V2[0],t0);
				/* V2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
				fp18_muln_fp3(t0,_S1[3],P_0);
				fp18_muln_fp3(t1,_P1[3],S_0);
				fp18_subc_low(t0,t0,t1);
				fp18_rdc(V2[1],t0);
				fp18_mul_dxs(V2[1],V2[1],beta);
				/* V2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
				fp18_muln_fp3(t0,_P1[4],S_0);
				fp18_muln_fp3(t1,_S1[4],P_0);
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
		fp3_free(alpha);	
		fp3_free(vt1);
		fp3_free(vt2);
		fp3_free(vt3);
		fp3_free(w2);
		fp3_free(w13);
		fp18_free(beta);
		fp18_free(gamma1);
		fp18_free(gamma2);
		fp18_free(delta);
		fp18_free(S_0);
		fp18_free(P_0);
		for(i=0;i<5;i++){
			fp3_free(_S1[i]);
			fp3_free(_P1[i]);
		}
		dv18_free(t0);
		dv18_free(t1);
}
}


void opt_en_tw_ri_kss(fp18_t F,ep3_t Q,ep_t P,fp3_t A,fp3_t B,bn_t m,int version){
	int i=0,sign=0;
	fp18_t initial_data[7],V2[3],xp,yp,w_n1,w_2n,w_21,w_11,k0,k1,x_q1,y_q1,x_mq,y_mq;
	fp3_t a,d,V1[8],_t0,_t1,_t2;
	bn_t l;

/*========================================================*/
/*                        Init                            */
/*========================================================*/
	TRY{
		bn_null(l);
		fp18_null(k0);
		fp18_null(k1);
		fp18_null(xp);
		fp18_null(yp);
		fp18_null(x_q1);
		fp18_null(y_q1);
		fp18_null(x_mq);
		fp18_null(y_mq);
		fp18_null(w_n1);
		fp18_null(w_2n);
		fp18_null(w_21);
		fp18_null(w_11);
		for(i=0;i<7;i++)	
			fp18_null(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_null(V2[i]);
		}
		fp3_null(a);
		fp3_null(d);
		fp3_null(_t0);
		fp3_null(_t1);
		fp3_null(_t2);
		for(i=0;i<8;i++){
			fp3_null(V1[i]);
		}
/*========================================================*/	
		bn_new(l);
		fp18_new(k0);
		fp18_new(k1);
		fp18_new(w_n1);
		fp18_new(w_2n);
		fp18_new(w_21);
		fp18_new(w_11);
		for(i=0;i<7;i++)	
			fp18_new(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_new(V2[i]);
		}
		fp3_new(a);
		fp3_new(d);
		fp18_new(xp);
		fp18_new(yp);
		fp18_new(x_q1);
		fp18_new(y_q1);
		fp18_new(x_mq);
		fp18_new(y_mq);
		fp3_new(_t0);
		fp3_new(_t1);
		fp3_new(_t2);
		for(i=0;i<8;i++){
			fp3_new(V1[i]);
		}
/*========================================================*/
		fp18_set_dig(F,1);
		if (ep_is_infty(P) || ep2_is_infty(Q)) return;
		fp3_copy(a,Q->x);//a=xq
		fp3_copy(d,Q->y);//d=yq
		
		fp18_zero(k0);
		fp_set_dig(k0[2][0][0],1);//w^2
		fp18_inv(k0,k0);
		fp18_zero(k1);
		fp_set_dig(k1[0][1][0],1);//w^3
		fp18_inv(k1,k1);

		fp18_zero(xp);
		fp18_zero(yp);
		fp_copy(xp[0][0][0],P->x);
		fp_copy(yp[0][0][0],P->y);
		
		fp18_mul(xp,xp,k0);
		fp18_mul(yp,yp,k1);

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
		fp3_set_dig(V1[3],1);//V1[3]=1
		fp3_dbl(V1[4],d);//the origin w(2,0)=2*yq
		/* the origin w(3,0)=3*xq^4+12*B*xq=3*xq*(xq^3+4*B) */	
		fp3_dbl(V1[5],a);
		fp3_add(V1[5],V1[5],a);//3*xq
		fp3_sqr(_t0,a);//t0=xq^2
		fp3_mul(_t0,_t0,a);//xq^3
		fp3_dbl(V1[2],B);//2B
		fp3_dbl(_t1,V1[2]);//4B
		fp3_add(V1[6],_t0,_t1);//xq^3+4B
		fp3_mul(V1[5],V1[5],V1[6]);
		/*the origin w(4,0)=4*yq*(xq^6+4B(5*xq^3-2*B))*/
		fp3_dbl(_t2,_t0);//2*xq^3
		fp3_add(_t2,_t2,_t2);
		fp3_add(_t2,_t2,_t0);//5*xq^3
		fp3_sub(V1[6],_t2,V1[2]);
		fp3_mul(V1[6],V1[6],_t1);
		fp3_sqr(_t0,_t0);//xq^6
		fp3_add(V1[6],V1[6],_t0);
		fp3_dbl(_t2,V1[4]);//4*yq
		fp3_mul(V1[6],V1[6],_t2);

		fp3_zero(V1[2]);//the orginal W(0,0)=0
		fp3_neg(V1[1],V1[3]);//the orginal W(-1,0)
		fp3_neg(V1[0],V1[4]);//the orginal W(-2,0)
		
	
	/* fill out the second vector of the block */
	
		fp18_set_dig(V2[0],1);// the original W(0,1)
		fp18_set_dig(V2[1],1);// the original W(1,1)
		/* the original W(2,1) =2*xq + xp - ((yp - yq)/(xp-xq))^2=2*a+xp-((yp-d)/(xp-a))^2*/
		/*fp6_zero(k0[1]);
		fp6_copy(k0[0],xp[0]);*/
		fp18_copy(k0,xp);

		fp_sub(k0[0][0][0],k0[0][0][0],a[0]);
		fp_sub(k0[0][2][0],k0[0][2][0],a[1]);
		fp_sub(k0[0][1][1],k0[0][1][1],a[2]);

		fp18_inv(k0,k0);//1/(xp-a)
		fp18_copy(k1,yp);
		fp_sub(k1[0][0][0],k1[0][0][0],d[0]);
		fp_sub(k1[0][2][0],k1[0][2][0],d[1]);
		fp_sub(k1[0][1][1],k1[0][1][1],d[2]);
		fp18_mul(k1,k1,k0);
		fp18_sqr(k1,k1);
		fp18_zero(k0);
		fp_dbl(k0[0][0][0],a[0]);
		fp_dbl(k0[0][2][0],a[1]);
		fp_dbl(k0[0][1][1],a[2]);
		fp18_add(k0,k0,xp);
		fp18_sub(V2[2],k0,k1);

		fp18_copy(w_n1,xp);
		fp18_neg(w_n1,w_n1);
		fp_add(w_n1[0][0][0],w_n1[0][0][0],a[0]);//w_n1=xq-xp=a-xp
		fp_add(w_n1[0][2][0],w_n1[0][2][0],a[1]);
		fp_add(w_n1[0][1][1],w_n1[0][1][1],a[2]);


		fp18_sqr(k1,w_n1);//k1=w_n1^2

		fp18_mul(k0,k1,k0);
		fp18_neg(w_2n,k0);

		fp18_copy(k0,yp);
		fp_add(k0[0][0][0],k0[0][0][0],d[0]);
		fp_add(k0[0][2][0],k0[0][2][0],d[1]);
		fp_add(k0[0][1][1],k0[0][1][1],d[2]);
		fp18_sqr(k0,k0);//(yp+d)^2
		fp18_add(w_2n,w_2n,k0);	
		fp18_copy(w_21,V2[2]);//the origin w(2,1)
		
		//V2[0]=V2[0]
		fp18_copy(V2[1],w_n1);//V2[1]=V2[1]*w_n1=1*w_n1=w_n1,new w(1,1)
		fp18_mul(V2[2],V2[2],k1);//new w(2,1)

		fp18_inv(k0,k1);//1/w_n1^2
		fp18_mul(w_2n,w_2n,k0);//w(-2,1)
		fp18_set_dig(w_n1,1);//w(1,-1)
		fp18_copy(w_11,V2[1]);//w(1,1)
		fp18_copy(w_21,V2[2]);//w(2,1)
		
    		fp3_inv(_t0,V1[4]);

		fp_copy(initial_data[0][0][0][0],_t0[0]);
		fp_copy(initial_data[0][0][2][0],_t0[1]);
		fp_copy(initial_data[0][0][1][1],_t0[2]);//w(2,0)^(-1)
		fp18_inv(initial_data[1],w_n1);//w_n1=1,w(-1,1)^(-1)
		fp18_inv(initial_data[2],w_2n);//w(-2,1)^(-1)
		fp18_inv(initial_data[3],w_21);//w(2,1)^(-1)
		fp18_inv(initial_data[4],w_11);//w(1,1)^(-1)
		fp3_sqr(_t0,V1[4]);
		fp3_mul(_t1,V1[3],V1[5]);
		fp18_zero(initial_data[5]);
		fp18_zero(initial_data[6]);
		
		fp_copy(initial_data[5][0][0][0],_t0[0]);
		fp_copy(initial_data[5][0][2][0],_t0[1]);
		fp_copy(initial_data[5][0][1][1],_t0[2]);
		fp_copy(initial_data[6][0][0][0],_t1[0]);
		fp_copy(initial_data[6][0][2][0],_t1[1]);
		fp_copy(initial_data[6][0][1][1],_t1[2]);

		switch(version){
			case 1:
				precomp(F,x_q1,y_q1,a,d,initial_data,V2,V1);
				DoubleAndAdd_2(initial_data,V2,V1,l);
				break;
			case 2:
				precomp_lz(F,x_q1,y_q1,a,d,initial_data,V2,V1);
				DoubleAndAdd_2_lz(initial_data,V2,V1,l);
				break;
			default:
				DoubleAndAdd_2(initial_data,V2,V1,l);
				break;
		}
	/*x_mq = xq - (V1[2] * V1[4]) / (V1[3] * V1[3])*/
		fp3_sqr(_t0,V1[3]);
		fp3_inv(_t1,_t0);
		fp3_mul(_t2,V1[2],V1[4]);	
		fp3_mul(_t1,_t1,_t2);
		fp3_sub(_t1,a,_t1);
		//notice that when m=r,V1[3]=0	
		fp18_zero(x_mq);
		fp_copy(x_mq[2][0][0],_t1[0]);
		fp_copy(x_mq[2][2][0],_t1[1]);
		fp_copy(x_mq[2][1][1],_t1[2]);


	/*y_mq=(V1[2]^2*V1[5]-V1[4]^2*V1[1])/4*yq*V1[3]^3*/
		fp3_mul(_t0,_t0,V1[3]);
		fp3_dbl(_t1,d);
		fp3_dbl(_t1,_t1);
		fp3_mul(_t0,_t0,_t1);
		fp3_inv(_t0,_t0);
		fp3_sqr(_t1,V1[2]);
		fp3_mul(_t1,_t1,V1[5]);
		fp3_sqr(_t2,V1[4]);
		fp3_mul(_t2,_t2,V1[1]);
		fp3_sub(_t1,_t1,_t2);
		fp3_mul(_t0,_t0,_t1);
		fp18_zero(y_mq);
		fp_copy(y_mq[0][1][0],_t0[0]);
		fp_copy(y_mq[0][0][1],_t0[1]);
		fp_copy(y_mq[0][2][1],_t0[2]);		

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
		fp_sub(x_q1[0][0][0],x_q1[0][0][0],P->x);
		fp18_mul(k0,k0,x_q1);
		fp18_neg(k0,k0);
		fp_sub(y_q1[0][0][0],y_q1[0][0][0],P->y);
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
		fp18_free(w_n1);
		fp18_free(w_2n);
		fp18_free(w_21);
		fp18_free(w_11);
		for(i=0;i<7;i++)	
			fp18_free(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_free(V2[i]);
		}
		fp3_free(a);
		fp3_free(d);
		fp18_free(xp);
		fp18_free(yp);
		fp3_free(_t0);
		fp3_free(_t1);
		fp3_free(_t2);
		fp18_free(x_q1);
		fp18_free(y_q1);
		fp18_free(x_mq);
		fp18_free(y_mq);
		for(i=0;i<8;i++){
			fp3_free(V1[i]);
		}
	}
}

void opt_en_tw_ri_kss_lz(fp18_t F,ep3_t Q,ep_t P,fp3_t A,fp3_t B,bn_t m,int version){
    	int i=0,sign=0;
	fp18_t initial_data[7],V2[3],xp,yp,w_n1,w_2n,w_21,w_11,k0,k1,x_q1,y_q1,x_mq,y_mq;
	fp3_t a,d,V1[8];
	dv3_t _t0,_t1;
    	fp3_t _t2,_t3,_t4;
	bn_t l;

/*========================================================*/
/*                        Init                            */
/*========================================================*/
	TRY{
		bn_null(l);
		fp18_null(k0);
		fp18_null(k1);
		fp18_null(xp);
		fp18_null(yp);
		fp18_null(x_q1);
		fp18_null(y_q1);
		fp18_null(x_mq);
		fp18_null(y_mq);
		fp18_null(w_n1);
		fp18_null(w_2n);
		fp18_null(w_21);
		fp18_null(w_11);
		for(i=0;i<7;i++)	
			fp18_null(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_null(V2[i]);
		}
		fp3_null(a);
		fp3_null(d);
		dv3_null(_t0);
		dv3_null(_t1);
		fp3_null(_t2);
		fp3_null(_t3);
		fp3_null(_t4);
		for(i=0;i<8;i++){
			fp3_null(V1[i]);
		}	
/*========================================================*/
		bn_new(l);
		fp18_new(k0);
		fp18_new(k1);
		fp18_new(w_n1);
		fp18_new(w_2n);
		fp18_new(w_21);
		fp18_new(w_11);
		for(i=0;i<7;i++)	
			fp18_new(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_new(V2[i]);
		}
		fp3_new(a);
		fp3_new(d);
		fp18_new(xp);
		fp18_new(yp);
		fp18_new(x_q1);
		fp18_new(y_q1);
		fp18_new(x_mq);
		fp18_new(y_mq);
		dv3_new(_t0);
		dv3_new(_t1);
		fp3_new(_t2);
		fp3_new(_t3);
		fp3_new(_t4);
		for(i=0;i<8;i++){
			fp3_new(V1[i]);
		}
/*========================================================*/	
		fp18_set_dig(F,1);
		if (ep_is_infty(P) || ep2_is_infty(Q)) return;
		fp3_copy(a,Q->x);//a=xq
		fp3_copy(d,Q->y);//d=yq
		
		fp18_zero(k0);
		fp_set_dig(k0[2][0][0],1);//w^2
		fp18_inv(k0,k0);
		fp18_zero(k1);
		fp_set_dig(k1[0][1][0],1);//w^3
		fp18_inv(k1,k1);

		fp18_zero(xp);
		fp18_zero(yp);
		fp_copy(xp[0][0][0],P->x);
		fp_copy(yp[0][0][0],P->y);
		
		fp18_mul(xp,xp,k0);
		fp18_mul(yp,yp,k1);

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

		fp3_set_dig(V1[3],1);//V1[3]=1
		fp3_dbl(V1[4],d);//the origin w(2,0)=2*yq
		/* the origin w(3,0)=3*xq^4+12*B*xq=3*xq*(xq^3+4*B) */
		fp3_dbl(V1[5],a);
		fp3_add(V1[5],V1[5],a);//3*xq
		fp3_sqr(_t2,a);//t0=xq^2
		fp3_mul(_t2,_t2,a);//xq^3
		fp3_dbl(V1[2],B);//2B
		fp3_dbl(_t3,V1[2]);//4B
		fp3_add(V1[6],_t2,_t3);//xq^3+4B
		fp3_mul(V1[5],V1[5],V1[6]);
		/*the origin w(4,0)=4*yq*(xq^6+4B(5*xq^3-2*B))*/
		
		fp3_dbl(V1[6],_t2);//2*xq^3
		fp3_dbl(V1[6],V1[6]);
		fp3_add(V1[6],V1[6],_t2);//5*xq^3
		fp3_sub(V1[6],V1[6],V1[2]);

		fp3_muln_low(_t0,V1[6],_t3);
		fp3_sqrn_low(_t1,_t2);//xq^6
		fp3_addc_low(_t0,_t1,_t0);
		fp3_rdc(V1[6],_t0);
		fp3_dbl(_t2,V1[4]);//4*yq
		fp3_mul(V1[6],V1[6],_t2);

		fp3_zero(V1[2]);//the orginal W(0,0)=0
		fp3_neg(V1[1],V1[3]);//the orginal W(-1,0)
		fp3_neg(V1[0],V1[4]);//the orginal W(-2,0)
	
	/* fill out the second vector of the block */
	
		fp18_set_dig(V2[0],1);// the original W(0,1)
		fp18_set_dig(V2[1],1);// the original W(1,1)
		/* the original W(2,1) =2*xq + xp - ((yp - yq)/(xp-xq))^2=2*a+xp-((yp-d)/(xp-a))^2*/
		/*fp6_zero(k0[1]);
		fp6_copy(k0[0],xp[0]);*/
		fp18_copy(k0,xp);

		fp_sub(k0[0][0][0],k0[0][0][0],a[0]);
		fp_sub(k0[0][2][0],k0[0][2][0],a[1]);
		fp_sub(k0[0][1][1],k0[0][1][1],a[2]);

		fp18_inv(k0,k0);//1/(xp-a)
		fp18_copy(k1,yp);
		fp_sub(k1[0][0][0],k1[0][0][0],d[0]);
		fp_sub(k1[0][2][0],k1[0][2][0],d[1]);
		fp_sub(k1[0][1][1],k1[0][1][1],d[2]);
		fp18_mul(k1,k1,k0);
		fp18_sqr(k1,k1);
		fp18_zero(k0);
		fp_dbl(k0[0][0][0],a[0]);
		fp_dbl(k0[0][2][0],a[1]);
		fp_dbl(k0[0][1][1],a[2]);
		fp18_add(k0,k0,xp);
		fp18_sub(V2[2],k0,k1);

		fp18_copy(w_n1,xp);
		fp18_neg(w_n1,w_n1);
		fp_add(w_n1[0][0][0],w_n1[0][0][0],a[0]);//w_n1=xq-xp=a-xp
		fp_add(w_n1[0][2][0],w_n1[0][2][0],a[1]);
		fp_add(w_n1[0][1][1],w_n1[0][1][1],a[2]);


		fp18_sqr(k1,w_n1);//k1=w_n1^2

		fp18_mul(k0,k1,k0);
		fp18_neg(w_2n,k0);

		fp18_copy(k0,yp);
		fp_add(k0[0][0][0],k0[0][0][0],d[0]);
		fp_add(k0[0][2][0],k0[0][2][0],d[1]);
		fp_add(k0[0][1][1],k0[0][1][1],d[2]);
		fp18_sqr(k0,k0);//(yp+d)^2
		fp18_add(w_2n,w_2n,k0);	
		fp18_copy(w_21,V2[2]);//the origin w(2,1)
		
		//V2[0]=V2[0]
		fp18_copy(V2[1],w_n1);//V2[1]=V2[1]*w_n1=1*w_n1=w_n1,new w(1,1)
		fp18_mul(V2[2],V2[2],k1);//new w(2,1)

		fp18_inv(k0,k1);//1/w_n1^2
		fp18_mul(w_2n,w_2n,k0);//w(-2,1)
		fp18_set_dig(w_n1,1);//w(1,-1)
		fp18_copy(w_11,V2[1]);//w(1,1)
		fp18_copy(w_21,V2[2]);//w(2,1)
		
		fp3_inv(_t2,V1[4]);

		fp_copy(initial_data[0][0][0][0],_t2[0]);
		fp_copy(initial_data[0][0][2][0],_t2[1]);
		fp_copy(initial_data[0][0][1][1],_t2[2]);//w(2,0)^(-1)
		fp18_inv(initial_data[1],w_n1);//w_n1=1,w(-1,1)^(-1)
		fp18_inv(initial_data[2],w_2n);//w(-2,1)^(-1)
		fp18_inv(initial_data[3],w_21);//w(2,1)^(-1)
		fp18_inv(initial_data[4],w_11);//w(1,1)^(-1)
		fp3_sqr(_t2,V1[4]);
		fp3_mul(_t3,V1[3],V1[5]);
		fp18_zero(initial_data[5]);
		fp18_zero(initial_data[6]);
		
		fp_copy(initial_data[5][0][0][0],_t2[0]);
		fp_copy(initial_data[5][0][2][0],_t2[1]);
		fp_copy(initial_data[5][0][1][1],_t2[2]);
		fp_copy(initial_data[6][0][0][0],_t3[0]);
		fp_copy(initial_data[6][0][2][0],_t3[1]);
		fp_copy(initial_data[6][0][1][1],_t3[2]);

		switch(version){
			case 1:
				precomp(F,x_q1,y_q1,a,d,initial_data,V2,V1);
				DoubleAndAdd_2(initial_data,V2,V1,l);
				break;
			case 2:
				precomp_lz(F,x_q1,y_q1,a,d,initial_data,V2,V1);
				DoubleAndAdd_2_lz(initial_data,V2,V1,l);
				break;
			default:
				DoubleAndAdd_2(initial_data,V2,V1,l);
				break;
		}
	/*x_mq = xq - (V1[2] * V1[4]) / (V1[3] * V1[3])*/
		fp3_sqr(_t2,V1[3]);
		fp3_inv(_t3,_t2);
		fp3_mul(_t4,V1[2],V1[4]);	
		fp3_mul(_t3,_t3,_t4);
		fp3_sub(_t3,a,_t3);
		//notice that when m=r,V1[3]=0	
		fp18_zero(x_mq);
		fp_copy(x_mq[2][0][0],_t3[0]);
		fp_copy(x_mq[2][2][0],_t3[1]);
		fp_copy(x_mq[2][1][1],_t3[2]);

	/*y_mq=(V1[2]^2*V1[5]-V1[4]^2*V1[1])/4*yq*V1[3]^3*/
		fp3_mul(_t2,_t2,V1[3]);
		fp3_dbl(_t3,d);
		fp3_dbl(_t3,_t3);
		fp3_mul(_t2,_t2,_t3);
		fp3_inv(_t2,_t2);
		fp3_sqr(_t3,V1[2]);
		fp3_mul(_t3,_t3,V1[5]);
		fp3_sqr(_t4,V1[4]);
		fp3_mul(_t4,_t4,V1[1]);
		fp3_sub(_t3,_t3,_t4);
		fp3_mul(_t3,_t3,_t2);
		fp18_zero(y_mq);
		fp_copy(y_mq[0][1][0],_t3[0]);
		fp_copy(y_mq[0][0][1],_t3[1]);
		fp_copy(y_mq[0][2][1],_t3[2]);		

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
		fp_sub(x_q1[0][0][0],x_q1[0][0][0],P->x);
		fp18_mul(k0,k0,x_q1);
		fp18_neg(k0,k0);
		fp_sub(y_q1[0][0][0],y_q1[0][0][0],P->y);
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
		fp18_free(w_n1);
		fp18_free(w_2n);
		fp18_free(w_21);
		fp18_free(w_11);
		for(i=0;i<7;i++)	
			fp18_free(initial_data[i]);
		for(i=0;i<3;i++){
			fp18_free(V2[i]);
		}
		fp3_free(a);
		fp3_free(d);
		fp18_free(xp);
		fp18_free(yp);
		dv3_free(_t0);
		dv3_free(_t1);
		fp3_free(_t2);
		fp3_free(_t3);
		fp3_free(_t4);
		fp18_free(x_q1);
		fp18_free(y_q1);
		fp18_free(x_mq);
		fp18_free(y_mq);
		for(i=0;i<8;i++){
			fp3_free(V1[i]);
		}
	}
}


