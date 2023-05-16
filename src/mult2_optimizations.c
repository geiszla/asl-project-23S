#include <math.h>
#define dbl_prec 53
#define binSize 45
// the following functions are defined in basefunctions.cpp
void twoMultFMA(const double a,const double b,double *pi_res, double *e_res);
int exponent(double d);
void accumulate(double p, double e, double *b, int sh, int l);
void renormalize(double *x,double* r, int n, int k);



/*
Optimization 0: Precompute 
*/
void mult2_0(double* x, double* y,double*pi, int n, int m, int r){
    int bins = r*dbl_prec/binSize+2;
		double B[bins];
    // get exponents
    int exp_x[n];
	int exp_y[m];
    for(int i=0;i<n;i++){
        exp_x[i] = exponent(x[i]);
    }
        for(int j=0;j<m;j++){
        exp_y[j] = exponent(y[j]);
    }
	// get sum of first exponents
	int e = exp_x[0] + exp_x[0];
	

    double B_start[bins];
	// initialize each Bin with starting value
    B_start[0] = ldexp(1.5,e+dbl_prec-1-binSize); // since 1.5*2^(e-(i+1)b+p-1) == 1.5*2^(e-b+p-1)*2^(-b*i)
    B[0] = B_start[0];
    double const C1 = ldexp(1.0, -binSize); // multiply in each iteration last value by 2^-b instead of 2^(-b*i)
	for (int i=1; i<bins;i++){
        B_start[i] = B_start[i-1]*C1; // 1.5*2^(e-(i+1)b+p-1)
		B[i] = B_start[i];
	}
	int j,l,sh;
	double p, err;
    int max_terms = fmin(n-1, r);
	for(int i=0; i<= max_terms; i++){
		for(j=0;j<= fmin(m-1,r-1-i); j++){
			//p = two_prod(x[i], y[j], &err); this leads to more similar result as CAMPARY
			twoMultFMA(x[i], y[j],&p,&err);
			l = e- exp_x[i] - exp_y[j]; 
			sh = floor(l/binSize); // bin of the first pair
			l = l- sh*binSize; // number of leading bits
			accumulate(p,err,B,sh,l); // add to correct bins
		}
		j-=1;
		if (j < m-1){ 
			p = x[i]*y[j];
			l = e- exp_x[i] - exp_y[j]; 
			sh = floor(l/binSize); 
			l = l- sh*binSize; 
			accumulate(p,0.,B,sh,l);
		}
	}
	for (int i=0;  i<bins;i++){
		B[i] = B[i]-B_start[i]; // B_i - 1.5*2^(e-(i+1)b+p-1)
	}

	
	renormalize(B,pi,bins,r);
}
/*
Optimization 0: Precompute 
Optimization 1: replace complex with simpler functions: fmin(a,b) -> a<b?a:b, floor -> (int),a/b -> a*b^-1

*/
void mult2_1(double* x, double* y,double*pi, int n, int m, int r){
    int bins = r*dbl_prec/binSize+2;
    int b_inv = 1/binSize;
	double B[bins];
    // get exponents
    int exp_x[n];
		int exp_y[m];
    for(int i=0;i<n;i++){
        exp_x[i] = exponent(x[i]);
    }
        for(int j=0;j<m;j++){
        exp_y[j] = exponent(y[j]);
    }
	// get sum of first exponents
	int e = exp_x[0] + exp_x[0];
	

    double B_start[bins];
	// initialize each Bin with starting value
    B_start[0] = ldexp(1.5,e+dbl_prec-1-binSize); // since 1.5*2^(e-(i+1)b+p-1) == 1.5*2^(e-b+p-1)*2^(-b*i)
    B[0] = B_start[0];
    double const C1 = ldexp(1.0, -binSize); // multiply in each iteration last value by 2^-b instead of 2^(-b*i)
	for (int i=1; i<bins;i++){
        B_start[i] = B_start[i-1]*C1; // 1.5*2^(e-(i+1)b+p-1)
		B[i] = B_start[i];
	}
	int j,l,sh;
	double p, err;
    int max_terms_x = (n-1)<r?(n-1):r;
	for(int i=0; i<= max_terms_x; i++){ // 
        int max_terms_y = ((m-1)<(r-1-i)?(m-1):(r-1-i));
		for(j=0;j<= max_terms_y; j++){ //(m-1)<(r-1-i)?(m-1):(r-1-i)
			//p = two_prod(x[i], y[j], &err); this leads to more similar result as CAMPARY
			twoMultFMA(x[i], y[j],&p,&err);
			l = e- exp_x[i] - exp_y[j]; 
			sh = (int)(l*b_inv); // bin of the first pair
			l = l- sh*binSize; // number of leading bits
			accumulate(p,err,B,sh,l); // add to correct bins
		}
		j-=1;
		if (j < m-1){ 
			p = x[i]*y[j];
			l = e- exp_x[i] - exp_y[j]; 
			sh = (int)(l*b_inv); 
			l = l- sh*binSize; 
			accumulate(p,0.,B,sh,l);
		}
	}
	for (int i=0;  i<bins;i++){
		B[i] = B[i]-B_start[i]; // B_i - 1.5*2^(e-(i+1)b+p-1)
	}

	
	renormalize(B,pi,bins,r);
}