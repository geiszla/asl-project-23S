#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctime>

extern "C"
{
  #include "basefunctions.h"
}

#define dbl_prec 53

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

void testvecsumerr(void (*implementation)(double *, int, double *) = vecSumErr) {
	int n = 38; //length of fp expansion: max of 38, otherwise matrix can't be nonoverlapping anymore and result will be trivial?

	double* f = (double *)malloc(n*sizeof(double));
	fill_matrix(f,n);

	// print
	// printf("Input:\n");
	// for (int i=0; i< n; i++){
	// 	printf("%f\n", f[i]);
	// }

	double* g = (double *)malloc(n*sizeof(double));

	implementation(f,n, g); // apply function
	fast_VecSumErr(f,n); // compute correct results (inplace)
	
	// printf("Output Algo 8:\n");
	// printf("OUR\t\t CAMPARY\n");
	for (int i=0; i< n; i++){
		// printf("%f\t %f \n", f[i], g[i]);
		assert(f[i]==g[i]);
	}

	free(f);
	free(g);
}

void testvecsumerrbranch(void (*implementation)(double *, int, int, double *) = vecSumErrBranch) {
	int n = 38; //length of fp expansion: max of 38, otherwise matrix can't be nonoverlapping anymore and result will be trivial?

	double* e = (double *)malloc(n*sizeof(double));
	fill_matrix(e, n);

	double* h = (double *)malloc(n*sizeof(double));
	
	implementation(e, n, n, h); // apply function
	fast_VecSumErrBranch(e, n, n); // compute correct results (inplace)
	
	// compare
	// printf("Output Algo 7:\n");
	// printf("OUR\t\t CAMPARY\n");
	for (int i=0; i< n; i++){
		// printf("%f\t %f \n", e[i], h[i]);
		assert(e[i]==h[i]);
	}

	free(e);
	free(h);
}
