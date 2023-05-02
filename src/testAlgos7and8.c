#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include <stdlib.h>
#define dbl_prec 53
////////////////////////////////////////copied from basefunctions.c//////////////////////////
void twoSum(double a, double b, double *s_res, double *e_res){
    double s = a+b;
    double t = s-b; 
    double e = (a-t) + (b-(s-t));
    // write ouput
    *s_res = s; *e_res = e;
    return;
}

void twoMultFMA(double a, double b,double *pi_res, double *e_res){
    double pi = a*b;
    double e = fma(a,b,-pi);
    *pi_res = pi; *e_res = e;
    return;

}

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
    return;

}


///////////////////////////////// copied from CAMPARY package/////////////////////////////
double FPadd_rn(const double x, const double y){
  return x + y;
}

double fast_two_sum(const double a, const double b, double* err){
  double s = FPadd_rn(a, b);
  double z = FPadd_rn(s, -a);
  *err = FPadd_rn(b, -z);
  return s;
}

void fast_VecSumErrBranch(double *x, int sX, int sR){
	int ptr = 0, i = 1;
	double e = x[0];

  	while(i<sX && ptr<sR){
		x[ptr] = fast_two_sum(e, x[i], &e); i++;
		if(e == 0.) e = x[ptr]; else ptr++;
  	}
  	if(ptr<sR && e!=0.){ x[ptr] = e; ptr++; }
		for(i=ptr; i<sR; i++) x[i] = 0.;
}

void fast_VecSumErr(double *x, int sX){
	double e;
	x[0] = fast_two_sum(x[0], x[1], &e);
	for(int i=2; i<sX-1; i++) x[i-1] = fast_two_sum(e, x[i], &e);
	x[sX-2] = fast_two_sum(e, x[sX-1], &x[sX-1]);
}


/**
Implementation of VecSumErrBranch algorithm (Algorithm 7)

Input: e vector size n (S-nonoverlapping), output vector size m
Output: f vector size m

**/
	double* vecSumErrBranch(double* e, int n, int m){
	double* err = (double *)malloc(n*sizeof(double));
	double* f = (double *)malloc(m*sizeof(double));
	int j = 0;
	err[0] = e[0];
	for (int i = 0; i <= n-2; i++) {
		twoSum(err[i], e[i+1], &f[j], &err[i+1]);
		if (err[i+1] != 0) {
			if (j >= m - 1){ // enough output terms
				return f;
			}
		j++;
		} else {
			err[i+1] = f[j];
		}
	}
	if (err[n-1] != 0 && j < m) {
		f[j] = err[n-1];
	}
	return f;
}


/**
Implementation of VecSumErr algorithm (Algorithm 8)

Input: f vector size n
Output: g vector size n

Only correct if n>2!!

Info: Probably faster to do it inplace in f then creating new g, but I tried to do it as similar as the algorithm.
**/
double* vecSumErr(double* f, int n){
	int m = n-1;
	double* err = (double *)malloc(n*sizeof(double));
	double* g = (double *)malloc(n*sizeof(double));
	err[0] = f[0];
	
	for (int i=0; i<= m-1; i++){
		twoSum(err[i],f[i+1],&g[i],&err[i+1]);
	}
	g[m] = err[m];
	free(err);
	return g;
}

// helpers
double ulp(double x){
                                                                   
   int exponent;                                                                       
                                                                                                                                                   
   frexp(x, &exponent);                                                            
                                                                                
   return ldexp(1.0,exponent); // 1.0*2^exponent = ulp
}

int exponent(double d)
{
  int result;
  frexp(d,&result);
  return result;
}

double pseudorand(double max)
{   
    return (max / RAND_MAX) * rand();
}
void fill_matrix(double * A, int n) {
	// assert P-nonoverlapping i.e. e_a[i-1] - e_a[i] >= p which implies also S nonoverlapping
	double m;
	int exp = 1023; // -1022 to +1023 in double
    for(int i=0; i < n; i++) {
    	m  = pseudorand(1);
        A[i] = ldexp(m,exp);
        exp -= 53;
    }
}

// test
int main()    
{  
	int n = 38; //length of fp expansion: max of 38, otherwise matrix can't be nonoverlapping anymore and result will be trivial?
	printf("n= %d\n",n);
	srand((unsigned) time(0)); //init random
	double* f = (double *)malloc(n*sizeof(double));
	fill_matrix(f,n);
	
	double* e = (double *)malloc(n*sizeof(double));
	fill_matrix(e, n);
	
	
	// print
	printf("Input:\n");
	for (int i=0; i< n; i++){
		printf("%f\n", f[i]);
		
	}
	
	double* g = vecSumErr(f,n); // apply function
	fast_VecSumErr(f,n); // compute correct results (inplace)
	
	double* h = vecSumErrBranch(e, n, n); // apply function
	fast_VecSumErrBranch(e, n, n); // compute correct results (inplace)
	
	// compare
	printf("Output Algo 7:\n");
	printf("OUR\t\t CAMPARY\n");
	for (int i=0; i< n; i++){
		printf("%f\t %f \n", e[i], h[i]);
		assert(e[i]==h[i]);
	}
	
	printf("Output Algo 8:\n");
	printf("OUR\t\t CAMPARY\n");
	for (int i=0; i< n; i++){
		printf("%f\t %f \n", f[i], g[i]);
		assert(f[i]==g[i]);
	}
	
	printf("Tests were successfull\n");
	

	
	free(f);
	free(g);
	return 1;
}