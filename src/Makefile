INCLUDE =  -I/opt/homebrew/Cellar/boost/1.81.0_1/include -I/opt/homebrew/Cellar/googletest/1.13.0/include

LIBPATH = -L/opt/homebrew/Cellar/googletest/1.13.0/lib

test:
	g++ unit_tests.cpp $(INCLUDE) $(LIBPATH) -lgtest -lgtest_main  -pthread -std=c++14 -o test
	
clean:
	rm test
