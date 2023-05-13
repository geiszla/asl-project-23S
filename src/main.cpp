/**
*      _________   _____________________  ____  ______
*     / ____/   | / ___/_  __/ ____/ __ \/ __ \/ ____/
*    / /_  / /| | \__ \ / / / /   / / / / / / / __/
*   / __/ / ___ |___/ // / / /___/ /_/ / /_/ / /___
*  /_/   /_/  |_/____//_/  \____/\____/_____/_____/
*
*  http://www.acl.inf.ethz.ch/teaching/fastcode
*  How to Write Fast Numerical Code 263-2300 - ETH Zurich
*  Copyright (C) 2019 
*                   Tyler Smith        (smitht@inf.ethz.ch) 
*                   Alen Stojanov      (astojanov@inf.ethz.ch)
*                   Gagandeep Singh    (gsingh@inf.ethz.ch)
*                   Markus Pueschel    (pueschel@inf.ethz.ch)
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program. If not, see http://www.gnu.org/licenses/.
*/
#ifdef _WIN32
#include <windows.h> // Include if under windows
#endif

#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "main.h"
#include "quaternion.h"
#include "timing.cpp"
#include "utils.h"

#ifdef __x86_64__
#include "tsc_x86.h"
#endif

using namespace std;

#define EPS (1e-3)

void kernel_base(quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N]) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
            A[i][j] = mul(x[i], y[j]);
        }
    }
}

void   register_functions();

/* Global vars, used to keep track of student functions */
template<typename T_compute_return, class... T_compute_args>
vector<compute_function<T_compute_return, T_compute_args...>> userFuncs;

vector<string> funcNames;
vector<int> funcFlops;
int numFuncs = 0;

/*
* Registers a user function to be tested by the driver program. Registers a
* string description of the function as well
*/
template<typename T_compute_return, class... T_compute_args>
void add_function(compute_function<T_compute_return, T_compute_args...> f, string name, int flops) {
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    funcFlops.push_back(flops);
    numFuncs++;
}


/*
* Checks the given function for validity. If valid, then computes and
* reports and returns the number of cycles required per iteration
*/
template<typename T_compute_return, class... T_compute_args>
double perf_test(compute_function<T_compute_return, T_compute_args...> f) {
    alignas(32) quaternion_t x[N];
    alignas(32) quaternion_t y[N];
    alignas(32) quaternion_t A[N][N];
    rands((double*) x, 4*N, 1);
    rands((double*) y, 4*N, 1);
    rands((double*) A, 4*N, N);

    // return 
    return measure_runtime(f, x, y, A);
}

int main(int argc, char **argv) {
  cout << "Starting program. ";
  double perf;
  int i;

  register_functions();

  if (numFuncs == 0){
    cout << endl;
    cout << "No functions registered - nothing for driver to do" << endl;
    cout << "Register functions by calling register_func(f, name)" << endl;
    cout << "in register_funcs()" << endl;

    return 0;
  }
  cout << numFuncs << " functions registered." << endl;
   
  //Check validity of functions.
  alignas(32) quaternion_t x[N];
  alignas(32) quaternion_t y[N];
  alignas(32) quaternion_t A[N][N];
  alignas(32) quaternion_t A_base[N][N];
  rands((double*) x, 4*N, 1);
  rands((double*) y, 4*N, 1);
  
  kernel_base(x, y, A_base);

  for (i = 0; i < numFuncs; i++) {
    auto f = userFuncs<void, quaternion_t[N], quaternion_t[N], quaternion_t[N][N]>[i];
    f(x, y, A);
    
    double error = nrm_sqr_diff((double*) A, (double*) A_base, 4*N*N);
    if (error > EPS) {
      cout << "\033[1;31m" << "The result of the " << i+1 << "th function is not correct." << "\033[0m" << std::endl;
    }
    rands((double*) A, 4*N, N);
  }

  for (i = 0; i < numFuncs; i++) {
    perf = perf_test(userFuncs<void, quaternion_t[N], quaternion_t[N], quaternion_t[N][N]>[i]);
    cout << endl << "Running: " << funcNames[i] << endl;
    cout << perf << " cycles" << endl;
  }

  return 0;
}