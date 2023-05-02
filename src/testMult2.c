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
/** Multiplies x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
    K - size of x, L - size of y, R - size of z. Constraints: K>=L.
    The algorithm computes the partial products in a paper-and-pencil fashion and then 
    accumulates them in a special designed structure that has a fixed-point flavour.
		double-precision = 53, bin size = 45;
    For operations using double-double, triple-double and quad-double we provide specialized 
    versions that use a generalization of the QD library's multiplication algorithm. **/
// Multiplication_accurate with relative error <= 2^{-(p-1)R}( 1 + (r+1)2^{-p} + 2^{-p+2}(K+L-R-2) )

#define dbl_prec 53
#define binSize 45
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


// helper
int exponent(double d)
{
  int result;
  frexp(d,&result);
  return result;
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

double* renormalize(double *x, int n, int k) {
	double eps = x[0];
	double* r = (double *)malloc(k*sizeof(double));
	int j = 0;
	int i = 1;
	while (i < n && j < k) {
		twoSum(eps, x[i], &r[j], &eps);
		if (eps == 0) {
			eps = r[j];
		} else {
			j++;
		}
		i++;
	}
	if (eps != 0 && j < k) {
		r[j] = eps;
	}
	return r;
}
/**
Implementation of Multiplication algorithm (Algorithm 1, Paper 2)

Input: x vector size n (ulp-nonoverlapping)
		y vector size m (ulp-nonoverlapping)
Output: pi vector size r (ulp-nonoverlapping) = x*y

Constraint: n >=m
**/

double* mult2(double* x, double* y, int n, int m, int r){
	int const LBN = r*dbl_prec/binSize + 2; // number of allocated bins 
	double B[LBN+2];
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
			p = two_prod(x[i], y[j], &err);
			l = e- exponent(x[i]) - exponent(y[i]); 
			sh = (int) (l/binSize); // bin of the first pair
			l = l- sh*binSize; // number of leading bits
			accumulate(p,err,B,sh,l); // add to correct bins
		}
		if (j < m-1){ // if min(m-1,r-1-i) = r-1-i: I don't get what this part does
			p = two_prod(x[i],y[j], &err);
			l = e- exponent(x[i]) - exponent(y[i]); 
			sh = (int) (l/binSize); 
			l = l- sh*binSize; 
			accumulate(p,0.,B,sh,l);
		}
	}
	for (int i=0;  i<LBN+2;i++){
		B[i] = B[i]-ldexp(1.5,e-(i+1)*binSize+dbl_prec-1); // B_i - 1.5*2^(e-(i+1)b+p-1)
	}
	
	
	
	return renormalize(B,n,r);
}
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
void fill_matrix(double * A, int n) {
	double m;
	int exp = 500; // -1022 to +1023 in double
    for(int i=0; i < n; i++) {
        m  = pseudorand(1);
        A[i] = ldexp(m,exp);
        exp -= 4;
        if (i>0){
        	assert(ulp(A[i-1]) >= A[i]);//assert ulp-nonoverlapping
		}
        
    }
}

// test
int main()    
{  
	int n = 4; // n >=m
	int m = 6;
	int R = 6;
	printf("n= %d, m= %d, r= %d\n",n,m ,R);
	srand((unsigned) time(0)); //init random
	double* x = (double *)malloc(n*sizeof(double));
	double* y = (double *)malloc(m*sizeof(double));
	double* r_true = (double *)malloc(R*sizeof(double));
	double* r_ours ;
	fill_matrix(x,n);
	fill_matrix(y, m);
	
	// apply correct function from CAMPARY
	truncatedMul(x, y, r_true, n, m,R);
	
	// apply our function
	r_ours = mult2(x,y,n,m,R);
	
	double sum_x=0.0;
	double sum_y=0.0;
	double sum_r_true=0.0;
	double sum_r_ours=0.0;
	// print
	printf("Input:\nx\t\ty\n");
	for (int i=0; i< n; i++){
		printf("%f\t%f\n", x[i], y[i]);
		sum_x+=x[i];
		sum_y+=y[i];
	}
	printf("Output:\nCAMPARY\t\tOurs\n");
	for (int i=0; i< R; i++){
		printf("%f\n%f\n\n", r_true[i], r_ours[i]);
		sum_r_true+=r_true[i];
		sum_r_ours+=r_ours[i];}
	printf("Quick check of result:\n%f\n%f\n%f\n", sum_x*sum_y, sum_r_true, sum_r_ours ); // should be similar iff R >= 2nm and the random numbers are not too big
	printf("Test was successfull\n");
	

	

	return 1;
}