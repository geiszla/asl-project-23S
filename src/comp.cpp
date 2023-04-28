#include <immintrin.h>
#include "common.h"
#include "complex.h"

// Function to print a 256-bit vector for debugging
void print_vector(__m256d v) {
    double* a = (double*)&v;
    printf("[%f, %f, %f, %f]\n", a[0], a[1], a[2], a[3]);
}

// Original function, unoptimized
void slow_performance1(quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N])
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      A[i][j] = mul(x[i], y[j]);
    }
  }
}

/*
 * Called by the driver to register your functions
 * Use add_function(func, description) to add your own functions
 */
void register_functions()
{
  add_function(&slow_performance1, "slow_performance1", 1);
}