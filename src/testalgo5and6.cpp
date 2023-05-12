#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"

#include <stdlib.h>
#include "basefunctions.cpp"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include <stdlib.h>

const double onedifference = pow(10,-16);
const double allonesindouble = 4.503599627370496 ;
///////////////////////////////// copied from CAMPARY package/////////////////////////////
double FPadd_rn(const double x, const double y){
  return x + y;
}
static inline double FPmul_rn(const double x, const double y){
  return x * y;
}
static inline double fma_d_rn_cpu(const double x, const double y, double xy) {
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

static inline double FPfma_rn(const double x, const double y, const double z){
	#ifdef FP_FAST_FMA
  //#warning cpu has fma
	  return fma(x, y, z);
	#else
   	return fma_d_rn_cpu(x, y, z);
	#endif
}
/* Computes fl(a*b) and err(a*b). */
static inline double two_prod(const double a, const double b, double &err){
	const double p = FPmul_rn(a, b);
 	err = FPfma_rn(a, b, -p);
 	return p;
}
static inline double two_sum(const double a, const double b, double &err){
  double s = FPadd_rn(a, b);
  double aa = FPadd_rn(s, -b);
  double bb = FPadd_rn(s, -aa);
  double da = FPadd_rn(a, -aa);
  double db = FPadd_rn(b, -bb);
  err = FPadd_rn(da, db);
  return s;
}
/* Computes fl(a+b) and err(a+b). Assumes |a| >= |b| */
static inline double fast_two_sum(const double a, const double b, double &err){
  double s = FPadd_rn(a, b);
  double z = FPadd_rn(s, -a);
  err = FPadd_rn(b, -z);
  return s;
}

/** Algorithm for normalizing the array x that contains random numbers. 
    After the first level	the result satisfies |x_i|<uls(x_{i-1}); S-nonoverlapping expansion.
		In the end, the result satisfies |x_i|<=ulp(x_{i-1}); ulp-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
// Renormalize_random

static inline void renorm_rand2L(int sX, int sR, double x[]){
	double pr;
	int i, j, ptr = 0;

  x[0] = two_sum(x[0], x[1], x[1]);
  for(i=2; i<sX; i++){
		pr = two_sum(x[i-1], x[i], x[i]);
    for(j=i-2; j>0; j--) pr = two_sum(x[j], pr, x[j+1]);
		x[0] = fast_two_sum(x[0], pr, x[1]);
	}

	i = 2;
  if(x[1] == 0.) pr = x[0];
  else { pr = x[1]; ptr++;}
  while(ptr<sR && i<sX){
    x[ptr] = fast_two_sum(pr, x[i], pr); i++;
    if(pr == 0.) pr = x[ptr]; else ptr++;
  }
  if(ptr<sR && pr!=0.){ x[ptr] = pr; ptr++; }
  for(i=ptr; i<sX; i++) x[i] = 0.;
}
static inline void merge(double const *x, double const *y, double *z,int K,int L){
  int i=0, j=0;
  for(int n=0; n<K+L; n++){
    if(i==K || (j<L && fabs(y[j])>fabs(x[i]))){ z[n] = y[j]; j++; }
    else{ z[n] = x[i]; i++; }
  }
}
static inline void certifiedAdd(const double *x, const double *y, double *r,int K,int L,int R){

    double f[K+L];
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


void testtwosum(){
		// test case modified as fast two sum gives sometimes results 
	for(double b = 1; b<100; b+=0.01){
		double a = randfrom(-onedifference*onedifference,+onedifference*onedifference );
		double c = randfrom(-onedifference*onedifference,+onedifference*onedifference );

		double err_ours,err_ref; double res;
		
		twoSum(a,b,&res,&err_ours);
		two_sum(a,b,err_ref);
		assert(err_ref==err_ours|| b ==0);
			
	}

}

void testfastfma(){
	// 2MultFMA
	for(double c = 1; c<100; c+=0.01){
		double a =  randfrom(-onedifference,+onedifference );
		double b =  randfrom(-onedifference,+onedifference);
		double res, err,error_ref;
		twoMultFMA(a,b,&res,&err);
		two_prod(a,b,error_ref);
		assert(error_ref==err||  b ==0);
	}
}

void testrenormalization(){
	// Test size 1 
	for(int c = 2; c<100; c+=5){
		double* renorm =  new double[c];
    double* solution =  new double[5];
    double*   sol_ref =  new double[5];
		for(int i = 0; i<c; i++){
			renorm[i] = allonesindouble;

		}
		renormalizationalgorithm(renorm,c,solution,1);
    renorm_rand2L(c,1,renorm);
    for(int n = 0; n<1; n++){
      double rt = renorm[n]; double st = solution[n];
      double a = 0;
      // test for reduction onto size 1 thus if failed thats the problem 
      assert(rt-st<0.00001);
    }
    double a=0;
	}

  for(int c = 2; c<20; c+=1){
		double* renorm =  new double[c];
    double* solution =  new double[5];
    double*   sol_ref =  new double[5];
		for(int i = 0; i<c; i++){
			renorm[i] = (allonesindouble*1024)/(8*i+1);

		}
		renormalizationalgorithm(renorm,c,solution,2);
    renorm_rand2L(c,2,renorm);
    for(int n = 0; n<2; n++){
      double rt = renorm[n]; double st = solution[n];
      double a = 0;
      // test for reduction onto size 2 thus if failed thats the problem 
      assert(rt-st<0.00001);
    }
    double a=0;
	}

  for(int c = 2; c<200; c+=1){
		double* renorm =  new double[c];
    double* solution =  new double[5];
    double*   sol_ref =  new double[5];
		for(int i = 0; i<c; i++){
			renorm[i] = (allonesindouble*1024)/(8*i+1);

		}
		renormalizationalgorithm(renorm,c,solution,3);
    renorm_rand2L(c,3,renorm);
    for(int n = 0; n<3; n++){
      double rt = renorm[n]; double st = solution[n];
      double a = 0;
      
      // test for reduction onto size 2 thus if failed thats the problem 
      assert(rt-st<0.00001);
    }
    double a=0;
	}
}


void testaddition(){
  for(int c = 2; c<200; c+=1){
		double* a =  new double[c];
    double* b =  new double[c];
    double*   sol =  new double[c];
    double*   sol_ref =  new double[c];
		for(int i = 0; i<c; i++){
			a[i] = 1;
      b[i] =1;
		}
		addition(a,b,sol,c,c,c);
    certifiedAdd(a,b,sol_ref,c,c,c);
    for(int i =0; i<c; i++){
      double sol1 = sol[i]; double sol2 = sol_ref[i];
      assert(sol1-(c*2 )<0.001);
    }
	}


}




int main()    
{  
  // testing framework 
  // two sum 

	testtwosum();
	testfastfma();
	testrenormalization();
  testaddition(); 

	


	printf("successfully runned all test cases");
}
