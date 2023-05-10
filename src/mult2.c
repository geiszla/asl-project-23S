#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include <stdlib.h>

#define dbl_prec 53
#define binSize 45
// compile with g++ -std=c++11 ./main.c
// error free transforms 
/*
*/
double mantissa(double d)
{
  int result;
  double m = frexp(d,&result);
  return m;
}
// Sorry a function can't start with a 2 therefore I took the nearest solution of it
void twoSum(double a, double b, double *s_res, double *e_res){
    double s = a+b;
    double t = s-b; 
    double e = (a-t) + (b-(s-t));
    // write ouput
    *s_res = s; *e_res = e;
    
}

void twoMultFMA(double a, double b,double *pi_res, double *e_res){
    double pi = a*b;
    double e = fma(a,b,-pi);
    *pi_res = pi; *e_res = e;
    

}


// call it with the array x and a array of the same size 
void vecSum(double *x, double *e_res){
    int length = sizeof(x)/sizeof(x[0]);
    assert(length == (sizeof(e_res)/sizeof(e_res[0])));
    double s[length]; 
    s[length-1] = x[length-1];
    for(int i = length-2; i>=0; i--){
        double s_tmp,e_tmp; 
        twoSum(x[i],s[i+1],&s_tmp, &e_tmp);
        s[i] = s_tmp; e_res[i] = e_tmp;
    }
    e_res[0]=s[0];

}

void deposit(double *p, double *b1, double *b2) {
	twoSum(*b1, *p, b1, p);
	*b2 = *b2 + *p;
}
/**
Implementation of Accumulate algorithm (Algorithm 2, Paper 2)

Input: 	double p, e
		vector<double> b
		int sh, l
Output: vector<double> b

**/

void accumulate(double p, double e, double *b, int sh, int l) {
	int c = dbl_prec - binSize - 1;
	if (l < binSize - 2*c - 1){
		deposit(&p, &b[sh], &b[sh+1]);
		deposit(&e, &b[sh+1], &b[sh+2]);
	} else if (l < binSize - c){
		deposit(&p, &b[sh], &b[sh+1]);
		twoSum(b[sh+1], e, &b[sh+1], &e);
		deposit(&e, &b[sh+2], &b[sh+3]);
	} else {
		twoSum(b[sh], p, &b[sh], &p);
		deposit(&p, &b[sh+1], &b[sh+2]);
		deposit(&e, &b[sh+2], &b[sh+3]);
	}
}

/**
Implementation of Renormalize algorithm (Algorithm 3, Paper 2)

Input: 	vector<double> x (size n)
		int n, k
Output: vector<double> r (ulp-nonoverlapping) (size k)
**/

void renormalize(double *x,double* r, int n, int k) {
	double eps = x[0];
	
	int j = 0;
	int i = 1;
	while (i < n && j < k) {
		twoSum(eps, x[i], &r[j], &eps);
		if (eps == 0) { // no overflow
			eps = r[j];
		} else {
			j++;
		}
		i++;
	}
	if (eps != 0 && j < k) {
		r[j] = eps;
	}
}

// helper
int exponent(double d)
{
  int result;
  frexp(d,&result);
  return result;
}
/**
Implementation of Multiplication algorithm (Algorithm 1, Paper 2)

Input: x vector size n (ulp-nonoverlapping)
		y vector size m (ulp-nonoverlapping)
Output: pi vector size r (ulp-nonoverlapping) = x*y

Constraint: n >=m
**/

double* mult2(double* x, double* y,double*pi, int n, int m, int r){
	int const LBN = r*dbl_prec/binSize; // number of allocated bins 
	double* B = (double *)malloc((LBN+2)*sizeof(double));
	// get sum of first exponents
	int e = exponent(x[0]) + exponent(y[0]);
	
	// initialize each Bin with starting value
	for (int i=0; i<LBN+2;i++){
		B[i] = ldexp(1.5,e-(i+1)*binSize+dbl_prec-1); // 1.5*2^(e-(i+1)b+p-1)
	}
	int j,l,sh;
	double p, err;
	for(int i=0; i<= fmin(n-1, r); i++){
		for(j=0;j<= fmin(m-1,r-1-i); j++){
			//p = two_prod(x[i], y[j], &err); this leads to more similar result as CAMPARY
			twoMultFMA(x[i], y[j],&p,&err);
			l = e- exponent(x[i]) - exponent(y[j]); 
			sh = floor(l/binSize); // bin of the first pair
			l = l- sh*binSize; // number of leading bits
			accumulate(p,err,B,sh,l); // add to correct bins
		}
		j-=1;
		if (j < m-1){ //  I don't get what this part does
			p = x[i]*y[j];
			l = e- exponent(x[i]) - exponent(y[j]); 
			sh = floor(l/binSize); 
			l = l- sh*binSize; 
			accumulate(p,0.,B,sh,l);
		}
	}
	for (int i=0;  i<LBN+2;i++){
		B[i] = B[i]-ldexp(1.5,e-(i+1)*binSize+dbl_prec-1); // B_i - 1.5*2^(e-(i+1)b+p-1)
	}

	
	//renorm_rand2L(B, pi, LBN+2,r);
	//fast_VecSumErrBranch(B,pi,LBN+2,r);
	renormalize(B,pi,LBN+2,r);
	return pi;
}




