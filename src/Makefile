ifeq ($(OS),Windows_NT)
	LIBPATH =
	INCLUDE =
	EXECUTABLE = test.exe
else
	INCLUDE =  -I/opt/homebrew/Cellar/boost/1.81.0_1/include -I/opt/homebrew/Cellar/googletest/1.13.0/include
	LIBPATH = -L/opt/homebrew/Cellar/googletest/1.13.0/lib
	EXECUTABLE = test
endif

test: unit_tests.cpp
	g++ unit_tests.cpp $(INCLUDE) $(LIBPATH) -lgtest -lgtest_main  -pthread -std=c++14 -o $(EXECUTABLE)
	./$(EXECUTABLE)
	
clean:
	rm test
