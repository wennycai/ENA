/*==========================================================================*/
/*===========                                              =================*/
/*                Original Elliptic Net Alg on BLS curves                    */
/*===========                                              =================*/
/*==========================================================================*/

#include"ena.h"

/*origin double and add alg*/
static void DoubleAndAdd(fp12_t initial_data[],fp12_t V2[],fp12_t V1[],bn_t m){
	fp12_t alpha,beta,gamma1,delta,S_0,P_0,t0;
	fp12_t S1[6],P1[6];
	int nb,i;

/*========================================================*/
/*                          Init                          */
/*========================================================*/
TRY{
	fp12_null(alpha);
	fp12_null(beta);
	fp12_null(gamma1);
	fp12_null(delta);
	fp12_null(S_0);
	fp12_null(P_0);
	fp12_null(t0);

	for(i=0;i<6;i++){
		fp12_null(S1[i]);
		fp12_null(P1[i]);
	}
/*========================================================*/	
	fp12_new(alpha);
	fp12_new(beta);
	fp12_new(gamma1);
	fp12_new(delta);
	fp12_new(S_0);
	fp12_new(P_0);
	fp12_new(t0);

	for(i=0;i<6;i++){
		fp12_new(S1[i]);
		fp12_new(P1[i]);
	}

/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

	nb = bn_bits(m);
	nb = nb - 1;
	
	//initial data contains the precomputed inverses
	fp12_copy(alpha,initial_data[0]);		
	fp12_copy(beta,initial_data[1]);// inverse of W(-1,1)
	fp12_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
	fp12_copy(delta,initial_data[4]); // inverse of W(1,1)	
	for(i=nb-1;i>=0;i--){
		int b=bn_get_bit(m,i);	
		fp12_sqr(S_0,V2[1]);//1S_12
		fp12_mul(P_0,V2[0],V2[2]);//1M_12
		for (int j = 0; j <= 5; j++)//6S_12,6M_12
		{
			fp12_sqr(S1[j],V1[j + 1]);
			fp12_mul(P1[j],V1[j],V1[j + 2]);
		}
		
		if (b == 0)
		{
			for (int j = 1; j <= 4; j++)//j=3 or 4
			{
			/* S[j - 1] * P1[j] - S[j] * P1[j - 1] */
				fp12_mul(V1[2*j-2],S1[j-1],P1[j]);
				fp12_mul(t0,S1[j],P1[j-1]);
				fp12_sub(V1[2*j-2],V1[2*j-2],t0);
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
				fp12_mul(V1[2*j-1],S1[j-1],P1[j+1]);	
				fp12_mul(t0,S1[j+1],P1[j-1]);
				fp12_sub(V1[2*j-1],V1[2*j-1],t0);
				fp12_mul_dxs(V1[2*j-1],V1[2*j-1],alpha);
			}

			/* V2[0] = (S1[1] * P_0 - P1[1] * S_0) * delta*/
			fp12_mul(V2[0],S1[1],P_0);
			fp12_mul(t0,P1[1],S_0);
			fp12_sub(V2[0],V2[0],t0);
			fp12_mul(V2[0],V2[0],delta);
			/*V2[1] = S1[2] * P_0 - P1[2] * S_0*/
			fp12_mul(V2[1],S1[2],P_0);
			fp12_mul(t0,P1[2],S_0);
			fp12_sub(V2[1],V2[1],t0);
			/* V2[2] = (S1[3] * P_0 - P1[3] * S_0) * beta */
			fp12_mul(V2[2],S1[3],P_0);
			fp12_mul(t0,P1[3],S_0);
			fp12_sub(V2[2],V2[2],t0);
			fp12_mul_dxs(V2[2],V2[2],beta);
		}
		else
		{
			for (int j = 1; j <= 4; j++)
			{
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
				fp12_mul(V1[2*j-2],S1[j-1],P1[j+1]);
				fp12_mul(t0,S1[j+1],P1[j-1]);
				fp12_sub(V1[2*j-2],V1[2*j-2],t0);
				fp12_mul_dxs(V1[2*j-2],V1[2*j-2],alpha);
			/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
				fp12_mul(V1[2*j-1],S1[j],P1[j+1]);		
				fp12_mul(t0,S1[j+1],P1[j]);
				fp12_sub(V1[2*j-1],V1[2*j-1],t0);
			}

			/* V2[0] = S1[2] * P_0 - P1[2] * S_0 */
			fp12_mul(V2[0],S1[2],P_0);
			fp12_mul(t0,P1[2],S_0);
			fp12_sub(V2[0],V2[0],t0);
			/* V2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
			fp12_mul(V2[1],S1[3],P_0);
			fp12_mul(t0,P1[3],S_0);
			fp12_sub(V2[1],V2[1],t0);
			fp12_mul_dxs(V2[1],V2[1],beta);
			/* V2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
			fp12_mul(V2[2],P1[4],S_0);
			fp12_mul(t0,S1[4],P_0);
			fp12_sub(V2[2],V2[2],t0);
			fp12_mul(V2[2],V2[2],gamma1);
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

	fp12_free(alpha);	
	fp12_free(beta);
	fp12_free(gamma1);
	fp12_free(delta);
	fp12_free(S_0);
	fp12_free(P_0);
	fp12_free(t0);

	for(i=0;i<6;i++){
		fp12_free(S1[i]);
		fp12_free(P1[i]);
	}
}
}

/*origin double and add alg with lazy reduction*/
static void DoubleAndAdd_lz(fp12_t initial_data[],fp12_t V2[],fp12_t V1[],bn_t m){
	fp12_t alpha,beta,gamma1,delta,S_0,P_0;
	fp12_t S1[6],P1[6];
	dv12_t t0,t1;
	int nb,i;

/*========================================================*/
/*                        Init                            */
/*========================================================*/
TRY{
	fp12_null(alpha);
	fp12_null(beta);
	fp12_null(gamma1);
	fp12_null(delta);
	fp12_null(S_0);
	fp12_null(P_0);
	dv12_null(t0);
	dv12_null(t1);

	for(i=0;i<6;i++){
		fp12_null(S1[i]);
		fp12_null(P1[i]);
	}
/*========================================================*/
	fp12_new(alpha);
	fp12_new(beta);
	fp12_new(gamma1);
	fp12_new(delta);
	fp12_new(S_0);
	fp12_new(P_0);
	dv12_new(t0);
	dv12_new(t1);

	for(i=0;i<6;i++){
		fp12_new(S1[i]);
		fp12_new(P1[i]);
	}

/*========================================================*/
/*                   Double-and-Add                       */
/*========================================================*/

	nb = bn_bits(m);
	nb = nb - 1;
	
	//initial data contains the precomputed inverses
	fp12_copy(alpha,initial_data[0]);		
	fp12_copy(beta,initial_data[1]);// inverse of W(-1,1)
	fp12_copy(gamma1,initial_data[2]);// inverse of W(-2,1)
	fp12_copy(delta,initial_data[4]); // inverse of W(1,1)	
	for(i=nb-1;i>=0;i--){
		int b=bn_get_bit(m,i);	
		fp12_sqr(S_0,V2[1]);//1S_12
		fp12_mul(P_0,V2[0],V2[2]);//1M_12
		for (int j = 0; j <= 5; j++)//6S_12,6M_12
		{
			fp12_sqr(S1[j],V1[j + 1]);
			fp12_mul(P1[j],V1[j],V1[j + 2]);
		}
		
		if (b == 0)
		{
			for (int j = 1; j <= 4; j++)//j=3 or 4
			{
			/* S[j - 1] * P1[j] - S[j] * P1[j - 1] */
				fp12_muln_low(t0,S1[j-1],P1[j]);
				fp12_muln_low(t1,S1[j],P1[j-1]);
				fp12_subc_low(t0,t0,t1);
				fp12_rdc(V1[2*j-2],t0);
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
				fp12_muln_low(t0,S1[j-1],P1[j+1]);	
				fp12_muln_low(t1,S1[j+1],P1[j-1]);
				fp12_subc_low(t0,t0,t1);
				fp12_rdc(V1[2*j-1],t0);
				fp12_mul_dxs(V1[2*j-1],V1[2*j-1],alpha);
			}

			/* V2[0] = (S1[1] * P_0 - P1[1] * S_0) * delta*/
			fp12_muln_low(t0,S1[1],P_0);
			fp12_muln_low(t1,P1[1],S_0);
			fp12_subc_low(t0,t0,t1);
			fp12_rdc(V2[0],t0);
			fp12_mul(V2[0],V2[0],delta);
			/*V2[1] = S1[2] * P_0 - P1[2] * S_0*/
			fp12_muln_low(t0,S1[2],P_0);
			fp12_muln_low(t1,P1[2],S_0);
			fp12_subc_low(t0,t0,t1);
			fp12_rdc(V2[1],t0);
			/* V2[2] = (S1[3] * P_0 - P1[3] * S_0) * beta */
			fp12_muln_low(t0,S1[3],P_0);
			fp12_muln_low(t1,P1[3],S_0);
			fp12_subc_low(t0,t0,t1);
			fp12_rdc(V2[2],t0);
			fp12_mul_dxs(V2[2],V2[2],beta);
		}
		else
		{
			for (int j = 1; j <= 4; j++)
			{
			/* (S[j - 1] * P1[j + 1] - S[j + 1] * P1[j - 1])*alpha */
				fp12_muln_low(t0,S1[j-1],P1[j+1]);
				fp12_muln_low(t1,S1[j+1],P1[j-1]);
				fp12_subc_low(t0,t0,t1);
				fp12_rdc(V1[2*j-2],t0);
				fp12_mul_dxs(V1[2*j-2],V1[2*j-2],alpha);
			/* S[j] * P1[j + 1] - S[j + 1] * P1[j] */
				fp12_muln_low(t0,S1[j],P1[j+1]);		
				fp12_muln_low(t1,S1[j+1],P1[j]);
				fp12_subc_low(t0,t0,t1);
				fp12_rdc(V1[2*j-1],t0);
			}

			/* V2[0] = S1[2] * P_0 - P1[2] * S_0 */
			fp12_muln_low(t0,S1[2],P_0);
			fp12_muln_low(t1,P1[2],S_0);
			fp12_subc_low(t0,t0,t1);
			fp12_rdc(V2[0],t0);
			/* V2[1] = (S1[3] * P_0 - P1[3] * S_0) * beta */
			fp12_muln_low(t0,S1[3],P_0);
			fp12_muln_low(t1,P1[3],S_0);
			fp12_subc_low(t0,t0,t1);
			fp12_rdc(V2[1],t0);
			fp12_mul_dxs(V2[1],V2[1],beta);
			/* V2[2] = (P1[4] * S_0 - S1[4] * P_0) * gamma1 */
			fp12_muln_low(t0,P1[4],S_0);
			fp12_muln_low(t1,S1[4],P_0);
			fp12_subc_low(t0,t0,t1);
			fp12_rdc(V2[2],t0);
			fp12_mul(V2[2],V2[2],gamma1);
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

	fp12_free(alpha);	
	fp12_free(beta);
	fp12_free(gamma1);
	fp12_free(delta);
	fp12_free(S_0);
	fp12_free(P_0);
	dv12_free(t0);
	dv12_free(t1);	

	for(i=0;i<6;i++){
		fp12_free(S1[i]);
		fp12_free(P1[i]);
	}
}
	
}

void opt_en_bls(fp12_t F,ep2_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version){
	int i=0,sign=0;
	fp12_t initial_data[7],V2[3],V1[8],xq,yq,w_n1,w_2n,w_21,w_11,k0,k1,B0,lambda;
	fp_t xp,yp;
	fp2_t a,d;	
	fp6_t _t0,_t1,_t2;
	bn_t l;

/*========================================================*/
/*                        Init                            */
/*========================================================*/
TRY{
	bn_null(l);
	fp12_null(k0);
	fp12_null(k1);	
	fp12_null(xq);
	fp12_null(yq);
	fp12_null(w_n1);
	fp12_null(w_2n);
	fp12_null(w_21);
	fp12_null(w_11);
	fp12_null(B0);
	for(i=0;i<7;i++)	
		fp12_null(initial_data[i]);
	for(i=0;i<3;i++){
		fp12_null(V2[i]);
	}
	fp2_null(a)
	fp2_null(d);
	fp_null(xp);
	fp_null(yp);
	fp12_null(lambda);
	fp6_null(_t0);
	fp6_null(_t1);
	fp6_null(_t2);
	for(i=0;i<8;i++){
		fp12_null(V1[i]);
	}	
/*========================================================*/	
	bn_new(l);
	fp12_new(k0);
	fp12_new(k1);	
	fp12_new(xq);
	fp12_new(yq);
	fp12_new(w_n1);
	fp12_new(w_2n);
	fp12_new(w_21);
	fp12_new(w_11);
	fp12_new(B0);
	for(i=0;i<7;i++)	
		fp12_new(initial_data[i]);
	for(i=0;i<3;i++){
		fp12_new(V2[i]);
	}
	fp2_new(a);
	fp2_new(d);
	fp_new(xp);
	fp_new(yp);
	fp12_new(lambda);
	fp6_new(_t0);
	fp6_new(_t1);
	fp6_new(_t2);
	for(i=0;i<8;i++){
		fp12_new(V1[i]);
	}
/*========================================================*/
	fp12_set_dig(F,1);
	if (ep_is_infty(P) || ep2_is_infty(Q)) return;
	fp2_copy(a,Q->x);//a=xq
	fp2_copy(d,Q->y);//d=yq
	fp_copy(xp,P->x);//xp
	fp_copy(yp,P->y);//yp

	fp12_zero(k0);
	fp_set_dig(k0[0][1][0],1);//w^2
	fp12_zero(k1);
	fp_set_dig(k1[1][1][0],1);//w^3
	fp12_inv(k0,k0);//1/w^2
	fp12_inv(k1,k1);//1/w^3
	fp12_zero(xq);
	fp12_zero(yq);
	fp2_copy(xq[0][0],a);
	fp2_copy(yq[0][0],d);
	fp12_mul(xq,xq,k0);
	fp12_mul(yq,yq,k1);	
	
	fp12_zero(B0);
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

	fp12_set_dig(V1[3],1);//V1[3]=1
	fp6_dbl(V1[4][1],yq[1]);//the origin w(2,0)=2*yq
	fp6_zero(V1[4][0]);
	/* the origin w(3,0)=3*xq^4+12*B*xq=3*xq*(xq^3+4*B) */
	fp6_zero(V1[5][1]);		
	fp6_dbl(V1[5][0],xq[0]);
	fp6_add(V1[5][0],V1[5][0],xq[0]);//3*xq
	fp12_sqr(k0,xq);//t0=xq^2
	fp12_mul(k0,k0,xq);//xq^3
	fp6_dbl(_t2,B0[0]);
	fp6_add(_t2,_t2,_t2);//4B
	fp6_zero(k1[1]);
	fp6_copy(k1[0],_t2);//k1=4B
	fp12_add(V1[6],k0,k1);
	fp12_mul(V1[5],V1[5],V1[6]);
	/*the origin w(4,0)=4*yq*(xq^6+4B(5*xq^3-2*B))*/
	fp6_zero(V1[6][1]);		
	fp6_dbl(_t1,k0[0]);//2*xq^3
	fp6_add(_t1,_t1,_t1);
	fp6_add(_t1,_t1,k0[0]);//5*xq^3
	fp6_dbl(V1[2][0],B0[0]);//2B
	fp6_sub(V1[6][0],_t1,V1[2][0]);
	fp12_mul(V1[6],V1[6],k1);
	fp12_sqr(k1,k0);//xq^6
	fp6_add(V1[6][0],V1[6][0],k1[0]);
	fp6_zero(V1[2][0]);
	fp6_dbl(V1[2][1],V1[4][1]);
	fp6_copy(_t2,V1[2][1]);//t2=4*yq
	fp12_mul(V1[6],V1[6],V1[2]);
	/*the orginal W(5,0)*/
	fp12_sqr(V1[7],V1[5]);
	fp12_mul(V1[7],V1[5],V1[7]);
	fp12_sqr(V1[2],V1[4]);
	fp12_mul(V1[2],V1[2],V1[4]);
	fp12_mul(V1[2],V1[2],V1[6]);
	fp12_sub(V1[7],V1[2],V1[7]);

	fp12_zero(V1[2]);//the orginal W(0,0)=0
	fp12_neg(V1[1],V1[3]);//the orginal W(-1,0)
	fp12_neg(V1[0],V1[4]);//the orginal W(-2,0)

/* fill out the second vector of the block */

	fp12_set_dig(V2[0],1);// the original W(0,1)
	fp12_set_dig(V2[1],1);// the original W(1,1)
	/* the original W(2,1) =2*xq + xp - ((yp - yq)/(xp-xq))^2*/
	fp6_neg(k0[0],xq[0]);
	fp_copy(k0[0][0][0],xp);
	fp12_inv(k0,k0);//1/(xp-xq)
	fp12_neg(k1,yq);
	fp_copy(k1[0][0][0],yp);
	fp12_mul(k1,k1,k0);
	fp12_sqr(k1,k1);
	fp6_dbl(k0[0],xq[0]);
	fp6_zero(k0[1]);
	fp_copy(k0[0][0][0],xp);//k0=2*xq + xp
	fp12_sub(V2[2],k0,k1);

	fp12_copy(w_n1,xq);
	fp_neg(w_n1[0][0][0],xp);//w_n1=xq-xp
	fp12_sqr(k1,w_n1);//k1=w_n1^2
	fp12_mul(k0,k1,k0);
	fp12_neg(w_2n,k0);
	fp12_copy(k0,yq);
	fp_copy(k0[0][0][0],yp);
	fp12_sqr(k0,k0);//(yp+y_q2)^2
	fp12_add(w_2n,w_2n,k0);	
	fp12_copy(w_21,V2[2]);//the origin w(2,1)
	
	//V2[0]=V2[0]
	fp12_copy(V2[1],w_n1);//V2[1]=V2[1]*w_n1=1*w_n1=w_n1,new w(1,1)
	fp12_mul(V2[2],V2[2],k1);//new w(2,1)

	fp12_inv(k0,k1);//1/w_n1^2
	fp12_mul(w_2n,w_2n,k0);//w(-2,1)
	fp12_set_dig(w_n1,1);//w(1,-1)
	fp12_copy(w_11,V2[1]);//w(1,1)
	
	fp12_inv(initial_data[0],V1[4]);//w(2,0)^(-1)
	fp12_inv(initial_data[1],w_n1);//w_n1=1,w(-1,1)^(-1)
	fp12_inv(initial_data[2],w_2n);//w(-2,1)^(-1)
	fp12_inv(initial_data[3],w_21);//w(2,1)^(-1)
	fp12_inv(initial_data[4],w_11);//w(1,1)^(-1)
	fp12_sqr(initial_data[5],V1[4]);
	fp12_mul(initial_data[6],V1[3],V1[5]);
				
	switch(version){
		case 1:
			DoubleAndAdd(initial_data,V2,V1,l);
			break;
		case 2:
			DoubleAndAdd_lz(initial_data,V2,V1,l);
			break;
		default:
			DoubleAndAdd(initial_data,V2,V1,l);
			break;
	}
//	fp12_inv(k0,V1[3]);
//	fp12_mul(F,k0,V2[1]);
	fp12_copy(F,V2[1]);
	/* f_{-a,Q}(P) = 1/f_{a,Q}(P). */
	if(sign==1){
		fp12_inv(F, F);
	}

	pp_exp_k12(F,F);

}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}	

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {

	bn_free(l);
	fp12_free(k0);
	fp12_free(k1);	
	fp12_free(xq);
	fp12_free(yq);
	fp12_free(w_n1);
	fp12_free(w_2n);
	fp12_free(w_21);
	fp12_free(w_11);
	fp12_free(B0);
	for(i=0;i<7;i++)	
		fp12_free(initial_data[i]);
	for(i=0;i<3;i++){
		fp12_free(V2[i]);
	}
	fp2_free(a);
	fp2_free(d);
	fp_free(xp);
	fp_free(yp);
	fp12_free(lambda);
	fp6_free(_t0);
	fp6_free(_t1);
	fp6_free(_t2);
	for(i=0;i<8;i++){
		fp12_free(V1[i]);
	}
}

}

void opt_en_bls_lz(fp12_t F,ep2_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version){
	int i=0,sign=0;
	fp12_t initial_data[7],V2[3],V1[8],xq,yq,w_n1,w_2n,w_21,w_11,k0,k1,B0,lambda;
	fp_t xp,yp;
	fp2_t a,d;	
	dv12_t _t0,_t1;
	bn_t l;

/*========================================================*/
/*                        Init                            */
/*========================================================*/
TRY{
	bn_null(l);
	fp12_null(k0);
	fp12_null(k1);	
	fp12_null(xq);
	fp12_null(yq);
	fp12_null(w_n1);
	fp12_null(w_2n);
	fp12_null(w_21);
	fp12_null(w_11);
	fp12_null(B0);
	for(i=0;i<7;i++)	
		fp12_null(initial_data[i]);
	for(i=0;i<3;i++){
		fp12_null(V2[i]);
	}
	fp2_null(a);
	fp2_null(d);
	fp_null(xp);
	fp_null(yp);
	fp12_null(lambda);
	dv12_null(_t0);
	dv12_null(_t1);
	for(i=0;i<8;i++){
		fp12_null(V1[i]);
	}	
/*========================================================*/
	bn_new(l);
	fp12_new(k0);
	fp12_new(k1);	
	fp12_new(xq);
	fp12_new(yq);
	fp12_new(w_n1);
	fp12_new(w_2n);
	fp12_new(w_21);
	fp12_new(w_11);
	fp12_new(B0);
	for(i=0;i<7;i++)	
		fp12_new(initial_data[i]);
	for(i=0;i<3;i++){
		fp12_new(V2[i]);
	}
	fp2_new(a);
	fp2_new(d);
	fp_new(xp);
	fp_new(yp);
	fp12_new(lambda);
	dv12_new(_t0);
	dv12_new(_t1);
	for(i=0;i<8;i++){
		fp12_new(V1[i]);
	}
/*========================================================*/
	fp12_set_dig(F,1);
	if (ep_is_infty(P) || ep2_is_infty(Q)) return;
	fp2_copy(a,Q->x);//a=xq
	fp2_copy(d,Q->y);//d=yq
	fp_copy(xp,P->x);//xp
	fp_copy(yp,P->y);//yp

	fp12_zero(k0);
	fp_set_dig(k0[0][1][0],1);//w^2
	fp12_zero(k1);
	fp_set_dig(k1[1][1][0],1);//w^3
	fp12_inv(k0,k0);//1/w^2
	fp12_inv(k1,k1);//1/w^3
	fp12_zero(xq);
	fp12_zero(yq);
	fp2_copy(xq[0][0],a);
	fp2_copy(yq[0][0],d);
	fp12_mul(xq,xq,k0);
	fp12_mul(yq,yq,k1);

	fp12_zero(B0);
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

	fp12_set_dig(V1[3],1);//V1[3]=1
	fp6_dbl(V1[4][1],yq[1]);//the origin w(2,0)=2*yq
	fp6_zero(V1[4][0]);
	/* the origin w(3,0)=3*xq^4+12*B*xq=3*xq*(xq^3+4*B) */
	fp12_sqr(k0,xq);//xq^2 mod p
	fp12_mul(k0,xq,k0);//xq^3 mod p
	fp12_dbl(k1,B0);//2B
	fp12_dbl(k1,k1);//4B
	fp12_add(V1[6],k0,k1);//xq^3+4B
	fp12_dbl(V1[5],xq);
	fp12_add(V1[5],V1[5],xq);
	fp12_mul(V1[5],V1[5],V1[6]);
	/*the origin w(4,0)=4*yq*(xq^6+20Bxq^3-8*B^2))*/
	fp12_sqrn_low(_t0,k0);//xq^6
	fp12_dbl(V1[6],k1);//8B
	fp12_muln_low(_t1,V1[6],B0);//8B^2
	fp12_subc_low(_t0,_t0,_t1);
	fp12_dbl(V1[6],V1[6]);//16B
	fp12_add(k1,k1,V1[6]);
	fp12_muln_low(_t1,k1,k0);
	fp12_addc_low(_t0,_t0,_t1);
	fp12_rdc(V1[6],_t0);
	fp12_dbl(V1[2],V1[4]);//4yp
	fp12_mul(V1[6],V1[6],V1[2]);
	/*the orginal W(5,0)*/
	fp12_sqr(V1[7],V1[5]);
	fp12_mul(V1[7],V1[5],V1[7]);
	fp12_sqr(V1[2],V1[4]);
	fp12_mul(V1[2],V1[2],V1[4]);
	fp12_mul(V1[2],V1[2],V1[6]);
	fp12_sub(V1[7],V1[2],V1[7]);
	fp12_zero(V1[2]);//the orginal W(0,0)=0
	fp12_neg(V1[1],V1[3]);//the orginal W(-1,0)
	fp12_neg(V1[0],V1[4]);//the orginal W(-2,0)

/* fill out the second vector of the block */

	fp12_set_dig(V2[0],1);// the original W(0,1)
	fp12_set_dig(V2[1],1);// the original W(1,1)
	/* the original W(2,1) =2*xq + xp - ((yp - yq)/(xp-xq))^2*/
	fp6_neg(k0[0],xq[0]);
	fp_copy(k0[0][0][0],xp);
	fp12_inv(k0,k0);//1/(xp-xq)
	fp12_neg(k1,yq);
	fp_copy(k1[0][0][0],yp);
	fp12_mul(k1,k1,k0);
	fp12_sqr(k1,k1);
	fp6_dbl(k0[0],xq[0]);
	fp6_zero(k0[1]);
	fp_copy(k0[0][0][0],xp);//k0=2*xq + xp
	fp12_sub(V2[2],k0,k1);

	fp12_copy(w_n1,xq);
	fp_neg(w_n1[0][0][0],xp);//w_n1=xq-xp
	fp12_sqr(k1,w_n1);//k1=w_n1^2
	fp12_mul(k0,k1,k0);
	fp12_neg(w_2n,k0);
	fp12_copy(k0,yq);
	fp_copy(k0[0][0][0],yp);
	fp12_sqr(k0,k0);//(yp+y_q2)^2
	fp12_add(w_2n,w_2n,k0);	
	fp12_copy(w_21,V2[2]);//the origin w(2,1)
	
	//V2[0]=V2[0]
	fp12_copy(V2[1],w_n1);//V2[1]=V2[1]*w_n1=1*w_n1=w_n1,new w(1,1)
	fp12_mul(V2[2],V2[2],k1);//new w(2,1)

	fp12_inv(k0,k1);//1/w_n1^2
	fp12_mul(w_2n,w_2n,k0);//w(-2,1)
	fp12_set_dig(w_n1,1);//w(1,-1)
	fp12_copy(w_11,V2[1]);//w(1,1)
	
	fp12_inv(initial_data[0],V1[4]);//w(2,0)^(-1)
	fp12_inv(initial_data[1],w_n1);//w_n1=1,w(-1,1)^(-1)
	fp12_inv(initial_data[2],w_2n);//w(-2,1)^(-1)
	fp12_inv(initial_data[3],w_21);//w(2,1)^(-1)
	fp12_inv(initial_data[4],w_11);//w(1,1)^(-1)
	fp12_sqr(initial_data[5],V1[4]);
	fp12_mul(initial_data[6],V1[3],V1[5]);
				
	switch(version){
		case 1:
			DoubleAndAdd(initial_data,V2,V1,l);
			break;
		case 2:
			DoubleAndAdd_lz(initial_data,V2,V1,l);
			break;
		default:
			DoubleAndAdd(initial_data,V2,V1,l);
			break;
	}
	fp12_inv(k0,V1[3]);
	fp12_mul(F,k0,V2[1]);
//	fp12_copy(F,V2[1]);
	/* f_{-a,Q}(P) = 1/f_{a,Q}(P). */
	if(sign==1){
		fp12_inv(F, F);
	}
	pp_exp_k12(F,F);

}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}	

/*========================================================*/
/*                        Free                            */
/*========================================================*/

	FINALLY {
	bn_free(l);
	fp12_free(k0);
	fp12_free(k1);	
	fp12_free(xq);
	fp12_free(yq);
	fp12_free(w_n1);
	fp12_free(w_2n);
	fp12_free(w_21);
	fp12_free(w_11);
	fp12_free(B0);
	for(i=0;i<7;i++)	
		fp12_free(initial_data[i]);
	for(i=0;i<3;i++){
		fp12_free(V2[i]);
	}
	fp2_free(a);
	fp2_free(d);
	fp_free(xp);
	fp_free(yp);
	fp12_free(lambda);
	dv12_free(_t0);
	dv12_free(_t1);
	for(i=0;i<8;i++){
		fp12_free(V1[i]);
	}
}
}
