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

// Eliminate all function calls
void eliminate_functions(quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N])
{
  for (int i = 0; i < 2; i++)
  {
    quaternion_t a = x[i];

    for (int j = 0; j < 2; j++)
    {
      quaternion_t b = y[j];
      
      double result_r, result_i, result_j, result_k;

      result_r = a.r * b.r - a.i * b.i - a.j * b.j - a.k * b.k;
      result_i = a.r * b.i + a.i * b.r + a.j * b.k - a.k * b.j;
      result_j = a.r * b.j - a.i * b.k + a.j * b.r + a.k * b.i;
      result_k = a.r * b.k + a.i * b.j - a.j * b.i + a.k * b.r;

      A[i][j] = { .r = result_r, .i = result_i, .j = result_j, .k = result_k };
    }
  }
}

// Vectorize code
void vectorization(quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N])
{
  for (int i = 0; i < 2; i++)
  {
    quaternion_t a = x[i];

    for (int j = 0; j < 2; j++)
    {
      // Load
      __m256d ar = _mm256_broadcast_sd(&a.r);
      __m256d b = _mm256_load_pd(&y[j].r);
      __m256d ai = _mm256_broadcast_sd(&a.i);
      __m256d aj = _mm256_broadcast_sd(&a.j);
      __m256d ak = _mm256_broadcast_sd(&a.k);
      __m256d mask1 = _mm256_set_pd(1, -1, 1, -1);
      __m256d mask2 = _mm256_set_pd(-1, 1, 1, -1);
      __m256d mask3 = _mm256_set_pd(-1, 1, -1, 1);

      // Shuffle
      __m256d b2 = _mm256_permute_pd(b, 0b0101);
      __m256d b3 = _mm256_permute4x64_pd(b, 0b01001110);
      __m256d b4 = _mm256_permute4x64_pd(b, 0b00011011);
      
      // Calculate
      __m256d x01, x02, x03, x04;

      x01 = _mm256_mul_pd(ar, b);
      x02 = _mm256_mul_pd(ai, b2);
      x03 = _mm256_mul_pd(aj, b3);
      x04 = _mm256_mul_pd(ak, b4);

      x02 = _mm256_mul_pd(x02, mask1);
      x01 = _mm256_add_pd(x01, x02);

      x04 = _mm256_mul_pd(x04, mask3);
      x03 = _mm256_add_pd(x03, x04);

      x03 = _mm256_mul_pd(x03, mask2);
      x01 = _mm256_add_pd(x01, x03);
      
      _mm256_store_pd(&A[i][j].r, x01);
    }
  }
}

// Unroll the inner loop
void unroll_inner_loop(quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N])
{
  for (int i = 0; i < 2; i++)
  {
    quaternion_t a = x[i];

    __m256d ar = _mm256_broadcast_sd(&a.r);
    __m256d ai = _mm256_broadcast_sd(&a.i);
    __m256d aj = _mm256_broadcast_sd(&a.j);
    __m256d ak = _mm256_broadcast_sd(&a.k);
    
    __m256d x01, x02, x03, x04, x05, x11, x12, x13, x14, x15;

    x02 = _mm256_load_pd(&y[0].r);
    x12 = _mm256_load_pd(&y[1].r);

    x01 = _mm256_mul_pd(ar, x02);
    x11 = _mm256_mul_pd(ar, x12);

    x04 = _mm256_permute_pd(x02, 0b0101);
    x03 = _mm256_mul_pd(ai, x04);
    x04 = _mm256_set_pd(1., -1., 1., -1.);
    x03 = _mm256_mul_pd(x03, x04);
    x01 = _mm256_add_pd(x01, x03);
    
    x14 = _mm256_permute_pd(x12, 0b0101);
    x13 = _mm256_mul_pd(ai, x14);
    x13 = _mm256_mul_pd(x13, x04);
    x11 = _mm256_add_pd(x11, x13);

    x04 = _mm256_permute4x64_pd(x02, 0b01001110);
    x03 = _mm256_mul_pd(aj, x04);
    
    x14 = _mm256_permute4x64_pd(x12, 0b01001110);
    x13 = _mm256_mul_pd(aj, x14);

    x05 = _mm256_permute4x64_pd(x02, 0b00011011);
    x04 = _mm256_mul_pd(ak, x05);
    x05 = _mm256_set_pd(-1., 1., -1., 1.);
    x04 = _mm256_mul_pd(x04, x05);
    x03 = _mm256_add_pd(x03, x04);

    x15 = _mm256_permute4x64_pd(x12, 0b00011011);
    x14 = _mm256_mul_pd(ak, x05);
    x14 = _mm256_mul_pd(x14, x15);
    x13 = _mm256_add_pd(x13, x14);
    
    x04 = _mm256_set_pd(-1., 1., 1., -1.);
    x03 = _mm256_mul_pd(x03, x04);
    x01 = _mm256_add_pd(x01, x03);

    x13 = _mm256_mul_pd(x13, x04);
    x11 = _mm256_add_pd(x11, x13);

    _mm256_store_pd(&A[i][0].r, x01);
    _mm256_store_pd(&A[i][1].r, x11);
  }
}

// Take out variables from the loop and rearrange operations
void miscellaneous_optimizations(quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N])
{
  for (int i = 0; i < 2; i++)
  {
    quaternion_t a = x[i];

    double a_r = x[i].r;
    double a_i = x[i].i;
    double a_j = x[i].j;
    double a_k = x[i].k;

    for (int j = 0; j < 2; j++)
    {
      // Load
      __m256d ar = _mm256_broadcast_sd(&a_r);
      __m256d b = _mm256_load_pd(&y[j].r);
      __m256d ai = _mm256_broadcast_sd(&a_i);
      __m256d aj = _mm256_broadcast_sd(&a_j);
      __m256d ak = _mm256_broadcast_sd(&a_k);

      // Shuffle
      __m256d b2 = _mm256_permute_pd(b, 0b0101);
      __m256d b3 = _mm256_permute4x64_pd(b, 0b01001110);
      __m256d b4 = _mm256_permute4x64_pd(b, 0b00011011);
      
      // Calculate
      __m256d x01, x02, x03, x04, x13, x14;

      x01 = _mm256_mul_pd(ar, b);
      x02 = _mm256_mul_pd(ai, b2);
      x01 = _mm256_addsub_pd(x01, x02);

      __m256d mask1 = _mm256_set_pd(-1, -1, 1, 1);
      x03 = _mm256_mul_pd(aj, b3);
      x03 = _mm256_mul_pd(x03, mask1);
      
      __m256d mask2 = _mm256_set_pd(1, -1, -1, 1);
      x04 = _mm256_mul_pd(ak, b4);
      x03 = _mm256_fmadd_pd(x04, mask2, x03);
      x01 = _mm256_addsub_pd(x01, x03);

      _mm256_store_pd(&A[i][j].r, x01);
      
      // double result_r, result_i, result_j, result_k;

      // result_r = (a.r * b.r - a.i * b.i) - (a.j * b.j + a.k * b.k);
      // result_i = (a.r * b.i + a.i * b.r) + (a.j * b.k - a.k * b.j);
      // result_j = (a.r * b.j - a.i * b.k) - (- a.j * b.r - a.k * b.i);
      // result_k = (a.r * b.k + a.i * b.j) + (- a.j * b.i + a.k * b.r);
  
      // A[i][j] = { .r = result_r, .i = result_i, .j = result_j, .k = result_k };
    }
  }
}

// The most optimized function (same code as in `miscellaneous_optimizations` above)
void maxperformance(quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N])
{
  for (int i = 0; i < 2; i++)
  {
    quaternion_t a = x[i];

    double a_r = x[i].r;
    double a_i = x[i].i;
    double a_j = x[i].j;
    double a_k = x[i].k;

    for (int j = 0; j < 2; j++)
    {
      // Load
      __m256d ar = _mm256_broadcast_sd(&a_r);
      __m256d b = _mm256_load_pd(&y[j].r);
      __m256d ai = _mm256_broadcast_sd(&a_i);
      __m256d aj = _mm256_broadcast_sd(&a_j);
      __m256d ak = _mm256_broadcast_sd(&a_k);

      // Shuffle
      __m256d b2 = _mm256_permute_pd(b, 0b0101);
      __m256d b3 = _mm256_permute4x64_pd(b, 0b01001110);
      __m256d b4 = _mm256_permute4x64_pd(b, 0b00011011);
      
      // Calculate
      __m256d x01, x02, x03, x04, x13, x14;

      x01 = _mm256_mul_pd(ar, b);
      x02 = _mm256_mul_pd(ai, b2);
      x01 = _mm256_addsub_pd(x01, x02);

      __m256d mask1 = _mm256_set_pd(-1, -1, 1, 1);
      x03 = _mm256_mul_pd(aj, b3);
      x03 = _mm256_mul_pd(x03, mask1);
      
      __m256d mask2 = _mm256_set_pd(1, -1, -1, 1);
      x04 = _mm256_mul_pd(ak, b4);
      x03 = _mm256_fmadd_pd(x04, mask2, x03);
      x01 = _mm256_addsub_pd(x01, x03);

      _mm256_store_pd(&A[i][j].r, x01);
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
  add_function(&eliminate_functions, "eliminate_functions", 1);
  add_function(&vectorization, "vectorization", 1);
  add_function(&unroll_inner_loop, "unroll_inner_loop", 1);
  add_function(&miscellaneous_optimizations, "miscellaneous_optimizations", 1);
  add_function(&maxperformance, "maxperformance", 1);
}