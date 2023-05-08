#include <iostream>
#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/number.hpp>
#include <gtest/gtest.h>
extern "C" {
#include "reference_solution.c"
}

using namespace boost::multiprecision;
typedef number<cpp_dec_float<1000> > big_float;

#define dbl_prec 53
#define binSize 45

double THRESHOLD = 0.9; // THIS NEEDS TO CHANGE AND DEPEND ON THE INPUT

double* create_ulp_non_overlaping_array(const int n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0, 1);
	double *x = (double*) malloc(sizeof(double) * n);
	for (int i = 0; i < n; i++) {
		x[i] = dist(gen);
		for (int j = n; j > i; j--) {
			x[i] *= 2;
		}
	}
	return x;
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

int main(int argc, char** argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
