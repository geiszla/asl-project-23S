#include <iostream>
#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/number.hpp>
#include <gtest/gtest.h>


#include "../lib/CAMPARY/Doubles/src_cpu/errorFreeTransf.h"
#include "../lib/CAMPARY/Doubles/src_cpu/renorm.h"
#include "../lib/CAMPARY/Doubles/src_cpu/addition.h"
#include "../lib/CAMPARY/Doubles/src_cpu/multiplication.h" 


extern "C" {
#include "basefunctions.c"
#include "mult2_optimizations.c"
}

// on windows compile with g++ unit_tests.cpp -lgtest -lgtest_main  -pthread -std=c++14 -o test.exe (or run "make test" to use makefile)
// for that you need to download the googleTest library and boost/multiprecision library. On windows i downloaded mingw-w64-i686-gtest and mingw-w64-boost from MSYS2 packages

using namespace std;

using namespace boost::multiprecision;
typedef number<cpp_dec_float<1000> > big_float;
#define dbl_prec 53
#define binSize 45



double get_ulp(double x){
	int exponent;
	frexp(x, &exponent);
	return ldexp(1.0,exponent); // 1.0*2^exponent = ulp
}


double* create_ulp_non_overlaping_array(const int n) {
	std::random_device rd;
	std::mt19937 gen(0); // set seed 
	std::uniform_real_distribution<double> dist(0.5, 1);
	double *x = (double*) malloc(sizeof(double) * n);
	double exponent = 200;
	for (int i = 0; i < n; i++) {
		x[i] = dist(gen)*pow(2,exponent);
		exponent-=53;
	}

	
	for (int i=1;i<n;i++) {
		EXPECT_LE(fabs(x[i]),get_ulp(x[i-1])); // assert ulp-nonoverlapping definition 2.6. (paper 2)
	}
	return x;
}
big_float get_error_tolerance(big_float x0, big_float y0, int N, int M, int R ){
	// error tolerance from the paper 2 (Proposition 3.5)
	return ldexp(fabs(x0*y0),-(dbl_prec-1)*R) * (1+ (R+1)*pow(2,-dbl_prec)+pow(2,-dbl_prec+1)*((-pow(2,-dbl_prec+1))/pow(1-pow(2,-dbl_prec+1),2)+(M+N-R-2)/(1-pow(2,-dbl_prec+1)))); 
}
big_float get_error_tolerance_certifiedMul(big_float r, int R ){
	// error tolerance from function description of certifiedMul in "/lib/CAMPARY/Doubles/src_cpu/multiplication.h"
	return pow(2,-(dbl_prec-1)*R)* (1+(fabs(r)+1)*pow(2,-dbl_prec)); 
}

big_float get_error_tolerance_truncatedMul(big_float r, int K, int L, int R ){
	// error tolerance from function description of truncatedMul in "/lib/CAMPARY/Doubles/src_cpu/multiplication.h"
	return pow(2,-(dbl_prec-1)*R)* (1+(fabs(r)+1)*pow(2,-dbl_prec) + pow(2,-dbl_prec+2)*(K+L-R-2)); 
}


big_float get_sum(double* x, const int n) {
	big_float sum = 0;
	for (int i = 0; i < n; i++) {
		sum += x[i];
	}
	return sum;
}
void print_array(double* x, int n) {

	for (int i = 0; i < n; i++) {
		cout << x[i] << "\n";}

}


TEST(ReferenceSolution, TruncatedMult) {
	const int K = 4;
	const int L = 4;
	const int R = 4;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	for (int i = 0; i < R; i++){
		r[i] = 0.0;
	}	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	truncatedMul<K,L,R>(x, y, r);
	

	big_float err_tolerance =  get_error_tolerance_truncatedMul(get_sum(r, R),K,L,R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));
	big_float rel_error = abs_error/get_sum(r, R);

	EXPECT_LT(rel_error, err_tolerance);
}
TEST(Multiplication, Mult2) {
	int K = 4;
	int L = 4;
	int R = 4;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	for (int i = 0; i < R; i++){
		r[i] = 0.0;
	}	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	mult2(x, y, r, K, L, R);
	

	big_float x0 = x[0];
	big_float y0 = y[0];
	
	big_float err_tolerance =  get_error_tolerance(x0, y0, K,L,R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));

/* 	//some printing TODO: remove this
	cout << "x:\n";
	print_array(x,K);
	cout << "y:\n";
	print_array(y,L);
	cout << "r:\n";
	print_array(r,R);
	cout << "abs error:" << abs_error << endl;
	cout << "err tolerance:" << err_tolerance << endl; */

	EXPECT_LT(abs_error, err_tolerance);
}
TEST(Multiplication, Mult2_error_tolerance_truncatedMul ) {
	int K = 4;
	int L = 4;
	int R = 4;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	for (int i = 0; i < R; i++){
		r[i] = 0.0;
	}	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	mult2(x, y, r, K, L, R);
	

	big_float err_tolerance =  get_error_tolerance_truncatedMul(get_sum(r, R),K,L,R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));
	big_float rel_error = abs_error/get_sum(r, R);

	EXPECT_LT(rel_error, err_tolerance);
}
TEST(Multiplication, Mult2_0 ) {
	int K = 4;
	int L = 4;
	int R = 4;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	for (int i = 0; i < R; i++){
		r[i] = 0.0;
	}	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	mult2_0(x, y, r, K, L, R);
	

	big_float err_tolerance =  get_error_tolerance_truncatedMul(get_sum(r, R),K,L,R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));
	big_float rel_error = abs_error/get_sum(r, R);

	EXPECT_LT(rel_error, err_tolerance);
}
TEST(Multiplication, Mult2_2 ) {
	int K = 4;
	int L = 4;
	int R = 4;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	for (int i = 0; i < R; i++){
		r[i] = 0.0;
	}	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	mult2_2(x, y, r, K, L, R);
	

	big_float err_tolerance =  get_error_tolerance_truncatedMul(get_sum(r, R),K,L,R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));
	big_float rel_error = abs_error/get_sum(r, R);

	EXPECT_LT(rel_error, err_tolerance);
}
TEST(Multiplication, Mult2_3 ) {
	int K = 4;
	int L = 4;
	int R = 4;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	for (int i = 0; i < R; i++){
		r[i] = 0.0;
	}	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	mult2_3(x, y, r, K, L, R);
	

	big_float err_tolerance =  get_error_tolerance_truncatedMul(get_sum(r, R),K,L,R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));
	big_float rel_error = abs_error/get_sum(r, R);

	EXPECT_LT(rel_error, err_tolerance);
}
TEST(ReferenceSolution, CertifiedMult) {
	const int K = 4;
	const int L = 4;
	const int R = 4;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	for (int i = 0; i < R; i++){
		r[i] = 0.0;
	}
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	certifiedMul<K,L,R>(x, y, r);
	

	
	big_float err_tolerance = get_error_tolerance_certifiedMul(get_sum(r, R),R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));
	big_float rel_error = abs_error/get_sum(r, R);

	EXPECT_LT(rel_error, err_tolerance);
}
int main(int argc, char** argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
