#include <cassert>
#include <cmath>


extern "C"
{
#ifdef __GNUC__
#include "reference.h"
#include "./basefunctions.c"
#include "./fourmult.c"

#else
#include "reference.h"
#include "basefunctions.h"
#include "./fourmult.c"
#endif
}

const double onedifference = pow(10,-16);
const double allonesindouble = 4.503599627370496 ;


/** Algorithm for normalizing the array x that contains random numbers. 
    After the first level	the result satisfies |x_i|<uls(x_{i-1}); S-nonoverlapping expansion.
		In the end, the result satisfies |x_i|<=ulp(x_{i-1}); ulp-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
// Renormalize_random


static inline void certifiedAdd(const double *x, const double *y, double *r,int K,int L,int R){

    double *f = new double[K+L];
    merge( x, y, f ,K,L);
    renorm_rand2L( K+L,R,f );
  
}





// from stackoverflow random generation
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}


void testtwosum(void (*implementation)(double, double, double *, double *) = twoSum){
		// test case modified as fast two sum gives sometimes results 
	for(double b = 1; b<100; b+=0.01){
		double a = randfrom(-onedifference*onedifference,+onedifference*onedifference );
		double c = randfrom(-onedifference*onedifference,+onedifference*onedifference );

		double err_ours,err_ref; double res;
		
		implementation(a,c,&res,&err_ours);
		two_sum(a,c,&err_ref);
		assert(err_ref==err_ours|| c ==0);
			
	}

}

void testfastfma(void (*implementation)(double, double, double *, double *) = twoMultFMA){
	// 2MultFMA
	for(double c = 1; c<100; c+=0.01){
		double a =  randfrom(-onedifference,+onedifference );
		double b =  randfrom(-onedifference,+onedifference);
		double res, err,error_ref;
		implementation(a,b,&res,&err);
		two_prod(a,b,&error_ref);
		assert(error_ref==err||  b ==0);
	}
}

void testrenormalization(void (*implementation)(double *, int, double *, int) = renormalizationalgorithm){
	// Test size 1 
	for(int c = 2; c<100; c+=5){
		double* renorm =  new double[c];
    double* solution =  new double[5];
    
		for(int i = 0; i<c; i++){
			renorm[i] = allonesindouble;

		}
		implementation(renorm,c,solution,1);
    renorm_rand2L(c,1,renorm);
    for(int n = 0; n<1; n++){
      double rt = renorm[n]; double st = solution[n];
     
      // test for reduction onto size 1 thus if failed thats the problem 
      assert(abs(rt-st)<0.00001);
    }
  
    delete[] renorm;
    delete[] solution;
	}

  for(int c = 2; c<20; c+=1){
		double* renorm =  new double[c];
    double* solution =  new double[5];
   
		for(int i = 0; i<c; i++){
			renorm[i] = (allonesindouble*1024)/(8*i+1);

		}
		implementation(renorm,c,solution,2);
    renorm_rand2L(c,2,renorm);
    
    for(int n = 0; n<2; n++){
      double rt = renorm[n];
      double st = solution[n];
     
      // test for reduction onto size 2 thus if failed thats the problem 
      assert(abs(rt-st)<0.00001);
    }
   
    delete[] renorm;
    delete[] solution;
	}

  for(int c = 3; c<200; c+=1){
		double* renorm =  new double[c];
    double* solution =  new double[5];
   
		for(int i = 0; i<c; i++){
			renorm[i] = (allonesindouble*1024)/(8*i+1);

		}
		implementation(renorm,c,solution,3);
  
    renorm_rand2L(c,3,renorm);
    for(int n = 0; n<3; n++){
      double rt = renorm[n];
      double st = solution[n];
      
      
      // test for reduction onto size 2 thus if failed thats the problem
      assert(abs(rt-st)<0.00001);
    }
   
    delete[] renorm;
    delete[] solution;
	}
}


void testaddition(void (*implementation)(double *, double *, double *, int, int, int) = addition){
  // assumption that the numbers are still non p-2 overlapping simmilar to the renormalization 
  // as addition is calling it and therefore it is kind of a natural implicit constraint
  for(int c = 2; c<20; c+=1){
		double* a =  new double[c];
    double* b =  new double[c];
    double*   sol =  new double[c]();
    double*   sol_ref =  new double[c]();
		for(int i = 0; i<c; i++){
			a[i] = (allonesindouble*1024)/(pow(8,2*i));
      b[i] = (allonesindouble*1024)/(pow(8,2*i+1));
		}
    // every 8 entries they go into the same "box"
    for(int i = 0; i<c; i+=8){
      for(int x = 0; x<8; x++){
        if(i+x < c){
          sol_ref[i/8] += (a[i+x]+b[i+x]);
        }
        
      }
    }
		implementation(a,b,sol,c,c,c);

    for(int i =0; i<c; i++){
      double ref = sol_ref[i]; double sol_our = sol[i];
      
      assert(abs(sol_our-ref) <0.001);
    }

    delete[] a;
    delete[] b;
    delete[] sol;
    delete[] sol_ref;
	}


   

}



void testmultiplication(void (*implementation)(double *, double *, double *, int, int, int) = multiplication){

  // reference implementation certifiedMul
  // our implementation implementation
  // test case 1
  for(int c = 3; c<10; c+=1){
    double* a =  new double[c];
    double* b =  new double[c];
    double*   sol =  new double[c]();
    double*   sol_ref =  new double[c]();
    for(int i = 0; i<c; i++){
      a[i] = (allonesindouble*1024)/(pow(8,2*i));
      b[i] = (allonesindouble*1024)/(pow(8,2*i+1));
    }
    renormalizationalgorithm(a, c, a, c);
    renormalizationalgorithm(b, c, b, c);
    implementation(a,b,sol,c,c,1);
    certifiedMul(c,c,1,a,b,sol_ref);
    for(int i =0; i<1; i++){
      double ref = sol_ref[i]; double sol_our = sol[i];
      assert(abs(sol_our-ref) <0.00001);
    } 
    delete[] a;
    delete[] b;
    delete[] sol;
    delete[] sol_ref;
  }

  // test case 2
  for(int c = 3; c<10; c+=1){
    double* a =  new double[c];
    double* b =  new double[c];
    double*   sol =  new double[c]();
    double*   sol_ref =  new double[c]();
    for(int i = 0; i<c; i++){
      a[i] = (allonesindouble*1024)/(pow(8,2*i));
      b[i] = (allonesindouble*1024)/(pow(8,2*i+1));
    }
    
    renormalizationalgorithm(a, c, a, c);
    renormalizationalgorithm(b, c, b, c);
    implementation(a,b,sol,c,c,2);
    certifiedMul(c,c,2,a,b,sol_ref);
    for(int i =0; i<2; i++){
      double ref = sol_ref[i]; double sol_our = sol[i];
      assert(abs(sol_our-ref) <0.00000001);
    } 
    delete[] a;
    delete[] b;
    delete[] sol;
    delete[] sol_ref;
  }
  // test case 3
  for(int c = 3; c<10; c+=1){
    double* a =  new double[c];
    double* b =  new double[c];
    double*   sol =  new double[c]();
    double*   sol_ref =  new double[c]();
    for(int i = 0; i<c; i++){
      a[i] = (allonesindouble*1024)/(pow(8,2*i));
      b[i] = (allonesindouble*1024)/(pow(8,2*i+1));
    }
    renormalizationalgorithm(a, c, a, c);
    renormalizationalgorithm(b, c, b, c);
    implementation(a,b,sol,c,c,3);
    certifiedMul(c,c,3,a,b,sol_ref);
    for(int i =0; i<3; i++){
      double ref = sol_ref[i]; double sol_our = sol[i];
      assert(abs(sol_our-ref) <0.00000001);
    } 
    delete[] a;
    delete[] b;
    delete[] sol;
    delete[] sol_ref;
  }


}


void testfourmultiplication(void (*implementation)(double *, double *, double *, double *, double *, double *,double *, double *, double *,double *, double *, double *,int, int, int) = fourtimesmultiplicationversion0){

 for(int c = 3; c<100; c+=1){
    double* a0 =  new double[c];
    double* b0 =  new double[c];
    double*   sol0 =  new double[c]();
    double*   sol_ref0 =  new double[c]();
    double* a1 =  new double[c];
    double* b1 =  new double[c];
    double*   sol1 =  new double[c]();
    double*   sol_ref1 =  new double[c]();
    double* a2 =  new double[c];
    double* b2 =  new double[c];
    double*   sol2 =  new double[c]();
    double*   sol_ref2 =  new double[c]();
    double* a3 =  new double[c];
    double* b3 =  new double[c];
    double*   sol3 =  new double[c]();
    double*   sol_ref3 =  new double[c]();

    for(int i = 0; i<c; i++){
      a0[i] = randfrom(-1.1,1.1);
      b0[i] = randfrom(-1.1,1.1);
      a1[i] = randfrom(-1.1,1.1);
      b1[i] = randfrom(-1.1,1.1);
      a2[i] = randfrom(-1.1,1.1);
      b2[i] = randfrom(-1.1,1.1);
      a3[i] = randfrom(-1.1,1.1);
      b3[i] = randfrom(-1.1,1.1);
    }
    renormalizationalgorithm(a0, c, a0, c);
    renormalizationalgorithm(b0, c, b0, c);
    renormalizationalgorithm(a1, c, a1, c);
    renormalizationalgorithm(b1, c, b1, c);
    renormalizationalgorithm(a2, c, a2, c);
    renormalizationalgorithm(b2, c, b2, c);
    renormalizationalgorithm(a3, c, a3, c);
    renormalizationalgorithm(b3, c, b3, c);
    implementation(a0,b0,a1,b1,a2,b2,a3,b3,sol0,sol1,sol2,sol3,c,c,4);
    certifiedMul(c,c,4,a0,b0,sol_ref0);
    certifiedMul(c,c,4,a1,b1,sol_ref1);
    certifiedMul(c,c,4,a2,b2,sol_ref2);
    certifiedMul(c,c,4,a3,b3,sol_ref3);
    for(int i =0; i<4; i++){
      double ref0 = sol_ref0[i]; double sol_our0 = sol0[i];
      double ref1 = sol_ref1[i]; double sol_our1 = sol1[i];
      double ref2 = sol_ref2[i]; double sol_our2 = sol2[i];
      double ref3 = sol_ref3[i]; double sol_our3 = sol3[i];
      assert(abs(sol_our0-ref0) <0.00000001);
      assert(abs(sol_our1-ref1) <0.00000001);
      assert(abs(sol_our2-ref2) <0.00000001);
      assert(abs(sol_our3-ref3) <0.00000001);
    }
    
    delete[] a0; 
    delete[] b0;
    delete[] sol0;
    delete[] sol_ref0;
    delete[] a1;
    delete[] b1;
    delete[] sol1;
    delete[] sol_ref1;
    delete[] a2;
    delete[] b2;
    delete[] sol2;
    delete[] sol_ref2;
    delete[] a3;
    delete[] b3;
    delete[] sol3;
    delete[] sol_ref3;

  }





}