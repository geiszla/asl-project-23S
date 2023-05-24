#include <cstdio>

#include "testalgo5and6.cpp"
#include "testAlgos7and8.cpp"

int main()
{
	// testing framework
	testmultiplication();
	testtwosum();
	testfastfma();
	//testvecsumerr();
	//testvecsumerrbranch();
	testrenormalization();
	testaddition();
	

	printf("successfully ran all test cases");
}
