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
