#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include <stdlib.h>
//#include "reference_solution.c"
#include "mult2.c"


// helper




// helpers
double ulp(double x){
                                                                   
   int exponent;                                                                       
                                                                                                                                                   
   frexp(x, &exponent);                                                            
                                                                                
   return ldexp(1.0,exponent); // 1.0*2^exponent = ulp
}
double pseudorand(double max)
{   
    return (max / RAND_MAX) * rand();
}

/**
	creates a matrix of size n with ulp-nonoverlapping fp-expansions
**/

double* create_matrix(int n) {
	double* A = (double *)malloc(n*sizeof(double));
	double m;
	int exp = 300; // -1022 to +1023 in double
    for(int i=0; i < n; i++) {
    	exp -= 53;
		
        m  = pseudorand(1); // random double 0-1
        A[i] = ldexp(m,exp);
	}
	double* B = (double *)malloc(n*sizeof(double));
  	//renorm_rand2L(A, B, n,n); // renormalize to ulp-nonoverlapping
  	fast_VecSumErrBranch(A,B,n,n);
  	return B;
	//A=renormalize(A,n,n);

}
// test
int main()    
{  
	int n = 16; // only m=n works
	int m = 16;
	int R = 16;
	printf("n= %d, m= %d, r= %d\n",n,m ,R);
	srand(1); //init seed
	double* x = create_matrix(n);
	double* y = create_matrix(m);
	double* r_true = (double *)malloc(R*sizeof(double));
	double* r_ours = (double *)malloc(R*sizeof(double));
	
	// apply correct function from CAMPARY
	certifiedMul(x, y, r_true, n, m,R);	

	double sum_x=0.0;
	double sum_y=0.0;
	double sum_r_true=0.0;
	double sum_r_ours=0.0;
	double dif=0.0;
	// print
	printf("Input:\nx\t\ty\n");
	for (int i=0; i< n; i++){
		printf("%lf*2^%d\t\t\t%lf*2^%d\n", mantissa(x[i]),exponent(x[i]), mantissa(y[i]), exponent(y[i]));
		sum_x+=x[i];
		sum_y+=y[i];
	}

	
	// apply our function
	mult2(x,y,r_ours,n,m,R);
	
	printf("Output:\nCAMPARY\t\t\t\tOurs\n");
	for (int i=0; i< R; i++){
		printf("%lf*2^%d\t\t\t%lf*2^%d\n", mantissa(r_true[i]),exponent(r_true[i]) , mantissa(r_ours[i]),exponent(r_ours[i]));
		sum_r_true+=r_true[i];
		sum_r_ours+=r_ours[i];
		dif += r_true[i] -r_ours[i];
		}
		

	double err_tolerance = ldexp(fabs(x[0]*y[0]),-(dbl_prec-1)*R) * (1+ (R+1)*pow(2,-dbl_prec)+pow(2,-dbl_prec+1)*((-pow(2,-dbl_prec+1))/pow(1-pow(2,-dbl_prec+1),2)+(m+n-R-2)/(1-pow(2,-dbl_prec+1)))); //error tolerance from the paper --> underflow to 0
	//err_tolerance = pow(2,-(dbl_prec-1)*R)*( 1 + (R+1)*pow(2,-dbl_prec) + pow(2,-dbl_prec+2)*(n+m-R-2) ); // error tolerance from the CAMPARY code
	printf("ErrorTolerance:\n%f\n", err_tolerance);
	printf("Quick check of result:\n%f(true)\n%f(CAMPARY)\n%f(Ours)\n", sum_x*sum_y, sum_r_true, sum_r_ours ); // should be similar iff R >= 2nm and the random numbers are not too big
	printf("Difference:\n%f\n",dif);
	printf("Test run successfull\n");
	

	

	return 1;
}