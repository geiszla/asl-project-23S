#include <math.h>


#define dbl_prec 53
#define binSize 45
//#include "FP_basicOp.h"
/*
 * FP_basicOp.h
 *
 * This file is part of CAMPARY Library
 *
 * Copyright (C) 2014 - 
 *
 * CAMPARY Library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * CAMPARY Library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MultiPrecGPU Library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, 
 * Boston, MA  02110-1301  USA
 */
 
/* Contributors: Valentina Popescu valentina.popescu@ens-lyon.fr
 *               Mioara Joldes joldes@laas.fr
 */


double FPadd_rn(const double x, const double y){
  return x + y;
}
static inline double FPmul_rn(const double x, const double y){
  return x * y;
}
static inline double FPfma_rn(const double x, const double y, double xy) {
  double C,px,qx,hx,py,qy,hy,lx,ly,fma;
  /*@ assert C == 0x1p27+1; */
	C = 0x1p27+1;

  px = x*C;
  qx = x-px;
  hx = px+qx;
  lx = x-hx;

  py = y*C;
  qy = y-py;
  hy = py+qy;
  ly = y-hy;

  fma = -x*y+hx*hy;
  fma += hx*ly;
  fma += hy*lx;
  fma += lx*ly;
  return fma;
}

/* Computes fl(a+b) and err(a+b). */
static inline double two_sum(const double a, const double b, double *err){
  double s = FPadd_rn(a, b);
  double aa = FPadd_rn(s, -b);
  double bb = FPadd_rn(s, -aa);
  double da = FPadd_rn(a, -aa);
  double db = FPadd_rn(b, -bb);
  *err = FPadd_rn(da, db);
  return s;
}
static inline double two_prod(const double a, const double b, double *err){
	const double p = FPmul_rn(a, b);
 	*err = FPfma_rn(a, b, -p);
 	return p;
}


double fast_two_sum(const double a, const double b, double* err){
  double s = FPadd_rn(a, b);
  double z = FPadd_rn(s, -a);
  *err = FPadd_rn(b, -z);
  return s;
}

void fast_VecSumErrBranch(double *x,double*r,  int sX, int sR){
	int ptr = 0, i = 1;
	double e = x[0];

  	while(i<sX && ptr<sR){
		r[ptr] = fast_two_sum(e, x[i], &e); i++;
		if(e == 0.) e = r[ptr]; else ptr++;
  	}
  	if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
		for(i=ptr; i<sR; i++) x[i] = 0.;
}

void fast_VecSumErr(double *x, int sX){
	double e;
	x[0] = fast_two_sum(x[0], x[1], &e);
	for(int i=2; i<sX-1; i++) x[i-1] = fast_two_sum(e, x[i], &e);
	x[sX-2] = fast_two_sum(e, x[sX-1], &x[sX-1]);
}

static inline void renorm_rand2L(const double *x, double *r, int sX, int sR){
  double pr, f[sX];
  int i, j, ptr = 0;

  f[0] = fast_two_sum(x[0], x[1], &f[1]);
  for(i=2; i<sX; i++){
    pr = fast_two_sum(f[i-1], x[i], &f[i]);
    for(j=i-2; j>0; j--) pr = fast_two_sum(f[j], pr, &f[j+1]);
    f[0] = fast_two_sum(f[0], pr, &f[1]);
  }

  i = 2;
  if(f[1] == 0.) pr = f[0];
  else { r[0] = f[0]; pr = f[1]; ptr++; }
  while(ptr<sR && i<sX){
    r[ptr] = fast_two_sum(pr, f[i], &pr); i++;
    if(pr == 0.) pr = r[ptr]; else ptr++;
  }
  if(ptr<sR && pr!=0.){ r[ptr] = pr; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

/** Multiplies x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
    K - size of x, L - size of y, R - size of z. Constraints: K>=L.
    The algorithm computes the partial products in a paper-and-pencil fashion and then 
    accumulates them in a special designed structure that has a fixed-point flavour.
		double-precision = 53, bin size = 45;
    For operations using double-double, triple-double and quad-double we provide specialized 
    versions that use a generalization of the QD library's multiplication algorithm. **/
// Multiplication_accurate with relative error <= 2^{-(p-1)R}( 1 + (r+1)2^{-p} + 2^{-p+2}(K+L-R-2) )


static inline void truncatedMul(const double *x, const double *y, double *r, int K, int L, int R){
	int const LBN = R*dbl_prec/binSize + 2;
  double B[LBN+2], lim[LBN+2];
  int i;

	int exp_x[(R+1<K)?R+1:K], exp_y[(R+1<L)?R+1:L];
	for(i=0; i<((R+1<K)?R+1:K); i++) frexp(x[i], &exp_x[i]);
	for(i=0; i<((R+1<L)?R+1:L); i++) frexp(y[i], &exp_y[i]);

	double factor = ldexp(1.0, -binSize); // 2^(-45)
	int exp_start = exp_x[0] + exp_y[0];
	lim[0] = ldexp(1.5, exp_start - binSize + dbl_prec-1); B[0] = lim[0];
	for(i=1; i<LBN+2; i++){ lim[i] = FPmul_rn(lim[i-1], factor); B[i] = lim[i]; }


		int j, l, sh;
		double p, e;
		for(i=0; i<(R<K?R:K); i++){
			for(j=0; j<(R-i<L?R-i:L); j++){
				l  = exp_start - (exp_x[i]+exp_y[j]);
				sh = (int)(l/binSize);
		  	l  = l - sh*binSize;
			  if(sh < LBN-1){
					p = two_prod(x[i], y[j], &e);
			    if(l < 30){ // binSize - 2*(dbl_prec-binSize-1) - 1){
						B[sh] = fast_two_sum(B[sh], p, &p);
						B[sh+1] = FPadd_rn(B[sh+1], p);

						B[sh+1] = fast_two_sum(B[sh+1], e, &e);
						B[sh+2] = FPadd_rn(B[sh+2], e);
					}else if(l < 37){ // binSize - (dbl_prec-binSize-1) - 1){
						B[sh] = fast_two_sum(B[sh], p, &p);
		        B[sh+1] = FPadd_rn(B[sh+1], p);

						B[sh+1] = fast_two_sum(B[sh+1], e, &e);
						B[sh+2] = fast_two_sum(B[sh+2], e, &e);
		        B[sh+3] = FPadd_rn(B[sh+3], e);
					}else{
						B[sh] = fast_two_sum(B[sh], p, &p);
						B[sh+1] = fast_two_sum(B[sh+1], p, &p);
		        B[sh+2] = FPadd_rn(B[sh+2], p);

						B[sh+2] = fast_two_sum(B[sh+2], e, &e);
		        B[sh+3] = FPadd_rn(B[sh+3], e);
		}}}}

	
	//computation of the error correction terms; using just simple multiplication
  if (R < L){
  	
  		double p;
  		int sh;
			for(i=0; i<=R; i++){
		  	sh = (int)((exp_start - (exp_x[i]+exp_y[R-i])) / binSize);
		    if(sh < LBN){
					p = FPmul_rn(x[i], y[R-i]);
					B[sh] = fast_two_sum(B[sh], p, &p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, &p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
		  }}
	}else if(R < K){
  		double p;
  		int sh;
		  for(i=0; i<L; i++){
				sh = (int)((exp_start - (exp_x[R-i]+exp_y[i])) / binSize);
		    if(sh < LBN){
		      p = FPmul_rn(x[R-i], y[i]);
					B[sh] = fast_two_sum(B[sh], p, &p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, &p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
		  }}
		
  }else{ 
  	if(R < K+L-1){
  		double p;
  		int sh;
		  for(i=R-L+1; i<K; i++){
				sh = (int)((exp_start - (exp_x[i]+exp_y[R-i])) / binSize);
		    if(sh < LBN){
		      p = FPmul_rn(x[i], y[R-i]);
					B[sh] = fast_two_sum(B[sh], p, &p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, &p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
			}}
	}
}
	/* unbias the B's */
	for (i=0; i<LBN; i++) B[i] = FPadd_rn(B[i], -lim[i]);
  fast_VecSumErrBranch(B, r, LBN,R);
}

/** Multiplies x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
    K - size of x, L - size of y, R - size of r.
    Computes all partial products based on scalar products and then accumulates 
    them in a special designed structure that has a fixed-point flavour.
    double-precision = 53, bin size = 45; **/
// Multiplication_accurate with relative error <= 2^{-(p-1)R}( 1 + (r+1)2^{-p} )

static inline void certifiedMul(const double *x, const double *y, double *r, int K, int L, int R){
	int const LBN = R*dbl_prec/binSize + 2;
	double B[LBN+2], lim[LBN+2];
	int i;
	
	int exp_x[K], exp_y[L];
	for(i=0; i<K; i++) frexp(x[i], &exp_x[i]);
	for(i=0; i<L; i++) frexp(y[i], &exp_y[i]);
	
	double factor = ldexp(1.0, -binSize); /* 2^(-45) */
	int exp_start = exp_x[0] + exp_y[0];
	lim[0] = ldexp(1.5, exp_start - binSize + dbl_prec-1); 
	
	B[0] = lim[0];
	for(i=1; i<LBN+2; i++){ 
		lim[i] = FPmul_rn(lim[i-1], factor); 
		B[i] = lim[i]; 
	}


	double p, e;
	int j, sh, l;
	for(i=0; i<K; i++){
	  for(j=0; j<L; j++){
	    l  = exp_start - (exp_x[i]+exp_y[j]);
	    sh = (int)(l/binSize);
	    l  = l - sh*binSize;
	    if(sh < LBN-1){
	      p = two_prod(x[i], y[j], &e);
	      if(l < 30){ // binSize - 2*(dbl_prec-binSize-1) - 1){
	        B[sh] = fast_two_sum(B[sh], p, &p);
	        B[sh+1] = FPadd_rn(B[sh+1], p);
	
	        B[sh+1] = fast_two_sum(B[sh+1], e, &e);
	        B[sh+2] = FPadd_rn(B[sh+2], e);
	      }else if(l < 37){ // binSize - (dbl_prec-binSize-1) - 1){
	        B[sh] = fast_two_sum(B[sh], p, &p);
	        B[sh+1] = FPadd_rn(B[sh+1], p);
	
	        B[sh+1] = fast_two_sum(B[sh+1], e, &e);
	        B[sh+2] = fast_two_sum(B[sh+2], e, &e);
	        B[sh+3] = FPadd_rn(B[sh+3], e);
	      }else{
	        B[sh] = fast_two_sum(B[sh], p, &p);
	        B[sh+1] = fast_two_sum(B[sh+1], p, &p);
	        B[sh+2] = FPadd_rn(B[sh+2], p);
	
	        B[sh+2] = fast_two_sum(B[sh+2], e, &e);
	        B[sh+3] = FPadd_rn(B[sh+3], e);
	}}}}
	

	/* unbias the B's */
	for (i=0; i<LBN; i++) B[i] = FPadd_rn(B[i], -lim[i]);
	fast_VecSumErrBranch(B, r, LBN,R);
}