#include <iostream>
#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/number.hpp>
#include <gtest/gtest.h>
extern "C" {
#include "reference_solution.c"
#include "mult2.c"
}



using namespace boost::multiprecision;
typedef number<cpp_dec_float<1000> > big_float;

#define dbl_prec 53
#define binSize 45

double THRESHOLD = 0.9; // THIS NEEDS TO CHANGE AND DEPEND ON THE INPUT


double get_ulp(double x){
	int exponent;
	frexp(x, &exponent);
	return ldexp(1.0,exponent); // 1.0*2^exponent = ulp
}

double* create_ulp_non_overlaping_array(const int n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0.5, 1);
	double *x = (double*) malloc(sizeof(double) * n);
	for (int i = 0; i < n; i++) {
		x[i] = dist(gen);
		for (int j = n; j > i; j--) {
			x[i] *= 2; //TODO: assert ulp-nonoverlapping
		}
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

big_float get_sum(double* x, const int n) {
	big_float sum = 0;
	for (int i = 0; i < n; i++) {
		sum += x[i];
	}
	return sum;
}


TEST(ReferenceSolution, TruncatedMult) {
	int K = 4;
	int L = 4;
	int R = 8;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	truncatedMul(x, y, r, K, L, R);
	
	big_float relative_error = get_sum(r, R) / (sum_x * sum_y);

	if (relative_error > 1) {
		relative_error = 1 / relative_error;
	}
	EXPECT_LT(THRESHOLD, relative_error);
}
TEST(Multiplication, Mult2) {
	int K = 4;
	int L = 4;
	int R = 8;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	r = mult2(x, y, r, K, L, R);
	

	big_float x0 = x[0];
	big_float y0 = y[0];
	
	big_float err_tolerance =  get_error_tolerance(x0, y0, K,L,R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));
	

	EXPECT_LT(abs_error, err_tolerance);
}
TEST(ReferenceSolution, CertifiedMult) {
	int K = 4;
	int L = 4;
	int R = 8;
	double *x = create_ulp_non_overlaping_array(K);
	double *y = create_ulp_non_overlaping_array(L);
	double *r = (double*) malloc(sizeof(double) * R);
	
	big_float sum_x = get_sum(x, K);
	big_float sum_y = get_sum(y, L);
	certifiedMul(x, y, r, K, L, R);
	
	
	big_float x0 = x[0];
	big_float y0 = y[0];
	
	big_float err_tolerance = get_error_tolerance(x0, y0, K,L,R);
	big_float abs_error = fabs((sum_x * sum_y)-get_sum(r, R));
	

	EXPECT_LT(abs_error, err_tolerance);
}
int main(int argc, char** argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
