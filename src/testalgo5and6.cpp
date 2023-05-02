#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include "reference_solution.c"

#include <stdlib.h>
#include "basefunctions.cpp"
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
        // test multiplication 
        double  a_arr[100];double b_arr[100];double r_arr[100];double r_arr_check[100];

        for(int i = 0; i<100; i++){
            a_arr[i]= pseudorand(1); 
            b_arr[i] =  pseudorand(1); 
        }
        multiplication(a_arr,b_arr,r_arr);
        //certifiedMul<100,100,100>(a_arr,b_arr,r_arr_check);

   }
   printf("all tests are successfull");

}
