#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include <stdlib.h>
const double trennung = 0.000000000000000001;// 10^-18
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
void multiplication(double *a, double *b, double *r,int sizea, int sizeb, int sizer){
    int k = sizea;
    double* err = (double *)malloc((sizea*sizea -1)*sizeof(double));
    double* pi_res = (double *)malloc(sizea*sizeof(double));
    double* r_ext = (double *)malloc(sizea*sizeof(double));
    twoMultFMA(a[0],b[0],&(r_ext[0]), &(err[0]));
    for(int n=1; n<sizea; n++){
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
        // write result into tmp
        vecSum(tmp2, tmp,1);
        // r_N = tmp[0]
        r_ext[n] = 0;
        // now write e 0:n^2 -1
        for(int b =0; b< n*n-1; b++){
            err[b]= tmp[b+1];
        }

        // now compute e[0:(n+1)^2 -1] <- e[0:n^2 + n -1],e[0:n]
        int count = 0;
        for(int b =0; b<=n^2 +n-1; b++){
            err[count] = err[b]; count++;
        }
        for(int b =0; b<=n; b++){
            err[count] = err[b]; count++;
        }
        // here compute line 9-11
        for ( int i =1; i<k; i++){
            r_ext[k] = r_ext[k] + a[i]*b[k-i];
        }
        // compute line 12 -14
         for ( int i =0; i<k*k; i++){
            r_ext[k] = r_ext[k] + err[i];
        }


        
        delete[] e_tmp;delete[] p; delete[] tmp;delete[] tmp2;
    }

    for(int i = 1; i<=sizea*sizea-1; i++){
        r_ext[k]= r_ext[k]+ err[i];
    }
    renormalizationalgorithm(r_ext,k+1,r,sizea);
    // return
}

// testing stuff

