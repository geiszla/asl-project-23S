#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include <stdlib.h>
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
void fast_VecSumErr(double *x, int sX){
	double e;
	x[0] = fast_two_sum(x[0], x[1], &e);
	for(int i=2; i<sX-1; i++) x[i-1] = fast_two_sum(e, x[i], &e);
	x[sX-2] = fast_two_sum(e, x[sX-1], &x[sX-1]);
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
double pseudorand(double max)
{   
    return (max / RAND_MAX) * rand();
}
void fill_matrix(double * A, int n) {
    for(int i=0; i < n; i++) {
        A[i] = pseudorand(1);
    }
}

// test
int main()    
{  
	int n = 5; //length of fp expansion
	printf("n= %d\n",n);
	srand((unsigned) time(0)); //init random
	double* f = (double *)malloc(n*sizeof(double));
	fill_matrix(f,n);
	
	
	// print
	printf("Input:\n");
	for (int i=0; i< n; i++){
		printf("%f\n", f[i]);
	}
	
	double* g = vecSumErr(f,n); // apply function
	fast_VecSumErr(f,n); // compute correct results (inplace)
	
	// compare
	printf("Output:\n");
	printf("OUR\t\t CAMPARY\n");
	for (int i=0; i< n; i++){
		printf("%f\t %f \n", f[i], g[i]);
		assert(f[i]==g[i]);
	}
	
	printf("Test was successfull\n");
	

	
	free(f);
	free(g);
	return 1;
}