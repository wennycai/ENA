#include <stdio.h>
#include "relic.h"
#include "relic_core.h"
#include "relic_fp_low.h"
#include "relic_fpx_low.h"
#include "relic_test.h"
#include "relic_bench.h"


//////////////////////////////////
void dv2_copy(dv2_t c,dv2_t a);
void fp2_rdc(fp2_t c, dv2_t a);
void fp3_rdc(fp3_t c, dv3_t a);
void fp12_rdc(fp12_t c, dv12_t a);
void fp18_rdc(fp18_t c, dv18_t a);
void fp6_muln_low(dv6_t c, fp6_t a, fp6_t b);
void fp12_muln_low(dv12_t c, fp12_t a, fp12_t b);
void fp12_sqrn_low(dv12_t c, fp12_t a);
void fp3_addc_low(dv3_t c,dv3_t a,dv3_t b);
void fp12_addc_low(dv12_t c,dv12_t a,dv12_t b);
void fp12_subc_low(dv12_t c,dv12_t a,dv12_t b);
void fp12_dbl(fp12_t c,fp12_t a);
void fp18_muln_low(dv18_t c, fp18_t a, fp18_t b);
void fp18_sqrn_low(dv18_t c, fp18_t a);
void fp18_addc_low(dv18_t c,dv18_t a,dv18_t b);
void fp18_subc_low(dv18_t c,dv18_t a,dv18_t b);
void fp18_dbl(fp18_t c,fp18_t a);


/*==========================================================================*/
/*               PAIRING VIA ELLIPTIC NET ON THE BLS CURVE                  */
/*==========================================================================*/

void opt_en_bls(fp12_t F,ep2_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);
void opt_en_bls_lz(fp12_t F,ep2_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);

void opt_en_cz_bls(fp12_t F,ep2_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);
void opt_en_cz_bls_lz(fp12_t F,ep2_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);

void opt_en_cz_ri_bls(fp12_t F,ep2_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);
void opt_en_cz_ri_bls_lz(fp12_t F,ep2_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);

void opt_en_tw_bls(fp12_t F,ep2_t Q,ep_t P,fp2_t A,fp2_t B,bn_t m,int version);
void opt_en_tw_bls_lz(fp12_t F,ep2_t Q,ep_t P,fp2_t A,fp2_t B,bn_t m,int version);

void opt_en_cz_tw_bls(fp12_t F,ep2_t Q,ep_t P,fp2_t A,fp2_t B,bn_t m,int version);
void opt_en_cz_tw_bls_lz(fp12_t F,ep2_t Q,ep_t P,fp2_t A,fp2_t B,bn_t m,int version);

void opt_en_tw_ri_bls(fp12_t F,ep2_t Q,ep_t P,fp2_t A,fp2_t B,bn_t m,int version);
void opt_en_tw_ri_bls_lz(fp12_t F,ep2_t Q,ep_t P,fp2_t A,fp2_t B,bn_t m,int version);


/*==========================================================================*/
/*               PAIRING VIA ELLIPTIC NET ON THE KSS CURVE                  */
/*==========================================================================*/

void opt_en_kss(fp18_t F,ep3_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);
void opt_en_kss_lz(fp18_t F,ep3_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);

void opt_en_cz_kss(fp18_t F,ep3_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);
void opt_en_cz_kss_lz(fp18_t F,ep3_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);

void opt_en_cz_ri_kss(fp18_t F,ep3_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);
void opt_en_cz_ri_kss_lz(fp18_t F,ep3_t Q,ep_t P,fp_t A,fp_t B,bn_t m,int version);

void opt_en_tw_kss(fp18_t F,ep3_t Q,ep_t P,fp3_t A,fp3_t B,bn_t m,int version);
void opt_en_tw_kss_lz(fp18_t F,ep3_t Q,ep_t P,fp3_t A,fp3_t B,bn_t m,int version);

void opt_en_cz_tw_kss(fp18_t F,ep3_t Q,ep_t P,fp3_t A,fp3_t B,bn_t m,int version);
void opt_en_cz_tw_kss_lz(fp18_t F,ep3_t Q,ep_t P,fp3_t A,fp3_t B,bn_t m,int version);

void opt_en_tw_ri_kss(fp18_t F,ep3_t Q,ep_t P,fp3_t A,fp3_t B,bn_t m,int version);
void opt_en_tw_ri_kss_lz(fp18_t F,ep3_t Q,ep_t P,fp3_t A,fp3_t B,bn_t m,int version);


/*==========================================================================*/
/*               SCALAR MULTIPLICATION VIA ELLIPTIC NET                     */
/*==========================================================================*/

void scalmul_net(ep_t mP,ep_t P,fp_t A,fp_t B,bn_t m,int version);
void scalmul_net_cz(ep_t mP,ep_t P,fp_t A,fp_t B,bn_t m,fp_t theta,int version);
