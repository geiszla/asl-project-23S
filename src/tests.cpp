#include <cstdio>

extern "C"
{
#include "reference.h"
#include "basefunctions.c"
#include "fourmult.c"
}

#include "testalgo5and6.cpp"
#include "testAlgos7and8.cpp"

int main()
{
	// testing framework
	testtwosum();
	testfastfma();
	testvecsumerr();
	testvecsumerrbranch();
	testrenormalization();
	testaddition();
	testmultiplication();
	testfourmultiplication();
	printf("successfully ran all test cases");
}
