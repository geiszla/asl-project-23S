#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctime>

extern "C"
{
  #include "reference.h"
  #include "basefunctions.h"
}

#define dbl_prec 53

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