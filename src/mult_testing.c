#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit/munit.h"
#include "mult2.c"
#include <stdlib.h>
#include <math.h>

static MunitResult test_mult2(const MunitParameter params[], void*data) {
	(void) params;
	(void) data;
	int n = 4; // only m=n works
	int m = 4;
	int R = 4;
	srand(1); //init seed
	double x[] = {1.0, 4.0, 16.0, 64.0};
	double y[] = {1.0, 16.0, 256.0, 4096.0};
	double* r_true = (double *)malloc(R*sizeof(double));
	double* r_ours = (double *)malloc(R*sizeof(double));
	
	// apply correct function from CAMPARY
	certifiedMul(x, y, r_true, n, m,R);	
	
	double sum_x=0.0;
	double sum_y=0.0;
	double sum_r_true=0.0;
	double sum_r_ours=0.0;
	double dif=0.0;
	mult2(x,y,r_ours,n,m,R);
	
	for (int i = 0; i < R; i++) {
		if (r_ours[i] != r_true[i]) {
			assert_double(r_ours[i], ==, r_true[i]);
			return MUNIT_FAIL;
		}
	}
	
	return MUNIT_OK;
}

static MunitTest test_suite_tests[] = {
	{
		(char*) "mult2",
		test_mult2,
		NULL,
		NULL,
		MUNIT_TEST_OPTION_NONE,
		NULL
	}
};

static const MunitSuite test_suite = {
	(char*) "",
	test_suite_tests,
	NULL,
	1,
	MUNIT_SUITE_OPTION_NONE
};


int main(int argc, char**argv) {
	return munit_suite_main(&test_suite, (void*) "µnit", argc, argv);
}