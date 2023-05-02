#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include "reference_solution.c"
#include <stdlib.h>

// compile with g++ -std=c++11 ./main.c
// error free transforms 
/*
*/

// Sorry a function can't start with a 2 therefore I took the nearest solution of it
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
    return;

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



/** implementation of Algorithm 6 renormalization
 **/

void renormalizationalgorithm(double *x, double *f, int m){
    int length = sizeof(x)/sizeof(x[0]);
    double e[length];
    vecSum(x,&e);
    f = vecSumErrBranch(e, length,m+1);
    for (int i =0; i<=m-2; i++){
        double *arr = vecSumErr(&f[i],m);
        for (int b=0; b<m; b++){f[b+i]=arr[b];}
    }

}


// testing stuff







// pseudorandom from stackoverflow
double pseudorand(double max)
{   
    srand((unsigned) time(0));
    return (max / RAND_MAX) * rand();
}


int main()    
{   
   for(int i = 0; i<10000; i++){
        double a = pseudorand(1); double b = pseudorand(1);double err_original,s_test,err_test; 
        two_sum(a,b, &err_original);
        twoSum(a,b, &s_test, &err_test);
        assert(err_original == err_test);

        two_prod(a,b, &err_original);
        twoMultFMA(a,b,&s_test, &err_test);
        assert(err_original == err_test);

   }
   printf("all tests are successfull");

}
