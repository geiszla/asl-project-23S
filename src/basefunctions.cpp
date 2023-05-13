#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include <stdlib.h>
const double trennung = 0.000000000000000001;// 10^-18

#define dbl_prec 53
#define binSize 45


// compile with g++ -std=c++11 ./main.cpp
// error free transforms 
/*
*/

// Sorry a function can't start with a 2 therefore I took the nearest solution of it
void twoSum(const double a,const  double b, double *s_res, double *e_res){
    double s = a+b;
    double t = s-b; 
    double e = (a-t) + (b-(s-t));
    // write ouput
    *s_res = s; *e_res = e;
    return;
}

void twoMultFMA(const double a,const double b,double *pi_res, double *e_res){
    double pi = a*b;
    double e = fma(a,b,-pi);
    *pi_res = pi; *e_res = e;
    return;

}


// call it with the array x and a array of the same size 
void vecSum(double *x, double *e_res, int in_out_size){
    int length = in_out_size;
    double *s = new double[length];
    s[length-1] = x[length-1];
    for(int i = length-2; i>=0; i--){
        double s_tmp,e_tmp; 
        twoSum(x[i],s[i+1],&s_tmp, &e_tmp);
        s[i] = s_tmp; e_res[i] = e_tmp;
    }
    e_res[0]=s[0];
    return;

}

/**
Implementation of VecSumErrBranch algorithm (Algorithm 7)

Input: e vector size n (S-nonoverlapping), output vector size m
Output: f vector size m

**/
double* vecSumErrBranch(double* e, int n, int m){
   double* err = new double[n];
   double* f = new double[n];
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
	double* err =  new double[n];
	double* g =  new double[n];
	err[0] = f[0];
	
	for (int i=0; i<= m-1; i++){
		twoSum(err[i],f[i+1],&g[i],&err[i+1]);
	}
	g[m] = err[m];
	
	return g;
}



/** implementation of Algorithm 6 renormalization
 **/

void renormalizationalgorithm(double x[],int size_of_x , double f[], int m){
    int length = size_of_x;
    double* e = new double[length];
    
    vecSum(x,e,size_of_x);
    double* f_tmp = vecSumErrBranch(e, length,m+1);
    for (int i =0; i<=m-2; i++){
        double *arr = vecSumErr(&f_tmp[i],m);
        for (int b=0; b<m; b++){
             double tmp=arr[b];
             f_tmp[b+i] = tmp;
            }
       double t= f_tmp[i];
       f[i] =t;
    }
     f[m-1]= f_tmp[m-1];
    delete[] e;
}

// camapry merge only for testing
static inline void merge(double const *x, double const *y, double *z,int K,int L){
  int i=0, j=0;
  for(int n=0; n<K+L; n++){
    if(i==K || (j<L && fabs(y[j])>fabs(x[i]))){ z[n] = y[j]; j++; }
    else{ z[n] = x[i]; i++; }
  }
}

/**Implementation of FP exansionaddition with k terms 
 * Input a and b of length k 
 * Output r of length k
*/
void addition(double *a, double *b, double *s, int length_a,int length_b, int length_result){
     double*  tmp = new double[length_a+length_b];
     merge(a,b,tmp,length_a,length_b);
     renormalizationalgorithm(tmp,length_a+length_b,s,length_result);
}

/**Implementation of FP exansion multiplication with k terms 
 * Input a and b of length k 
 * Output r of length k
*/
void multiplication(double *a, double *b, double *r){
    int length = sizeof(a)/sizeof(a[0]);int k = length;
    double* err = (double *)malloc((k*k+k-1)*sizeof(double));
    double* pi_res = (double *)malloc(length*sizeof(double));
    twoMultFMA(a[0],b[0],&(pi_res[0]), &(err[0]));
    for(int n=1; n<length; n++){
        double*  e_tmp = new double[n];
        double*  p = new double[n];
        for(int i = 0; i<=n; i++){
            twoMultFMA(a[i],b[n-i],&(p[i]), &(e_tmp[i]));
        }
        double*  tmp = new double[n*n+n];
        double*  tmp2 = new double[n*n+n+1];
        for(int b =0; b<=n; b++){
            tmp2[b]=p[b];
        }
        for(int b =n+1; b<n*n+n; b++){
            tmp2[b]=err[b-n];
        }
        vecSum(tmp2, tmp,1);
        r[n] = tmp[0];
        for(int b = 1; b<=n*n+n; b++){
            err[b-1] = tmp[b];
        }
        for(int b = 0; b<=n; b++){
            err[n*n+n+ b+1] = e_tmp[b];
        }
        
        delete[] e_tmp;delete[] p; delete[] tmp;delete[] tmp2;
    }

    for(int i = 1; i<=k-1; i++){
        r[k]= r[k]+ err[i];
    }
    renormalizationalgorithm(r,2,r,k);
    // return
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

