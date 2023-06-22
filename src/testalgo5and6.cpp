#include <stdio.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
/*
#include "vec_sum_optimizations.c"
#include "vec_sum_err_optimizations.c"
#include "vec_sum_err_branch_optimizations.c"
#include "renormalization_optimizations.c"
#include "addition_optimizations.c"
#include "multiplication_optimizations.c"
*/

const double onedifference = pow(10, -16);
const double allonesindouble = 4.503599627370496;

/** Algorithm for normalizing the array x that contains random numbers.
    After the first level	the result satisfies |x_i|<uls(x_{i-1}); S-nonoverlapping expansion.
    In the end, the result satisfies |x_i|<=ulp(x_{i-1}); ulp-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
// Renormalize_random

static inline void certifiedAdd(const double *x, const double *y, double *r, int K, int L, int R)
{

  double *f = new double[K + L];
  merge(x, y, f, K, L);
  renorm_rand2L(K + L, R, f);
}
void sort_double_array_descending(double *array, size_t size) {
  if (size <= 1) {
    return;
  }

  // Divide the array into two halves.
  size_t mid = size / 2;
  double *left_array = (double *) malloc(sizeof(double) * mid);
  double *right_array = (double *) malloc(sizeof(double) * (size - mid));

  for (size_t i = 0; i < mid; i++) {
    left_array[i] = array[i];
  }

  for (size_t i = mid; i < size; i++) {
    right_array[i - mid] = array[i];
  }

  // Sort the left and right halves recursively.
  sort_double_array_descending(left_array, mid);
  sort_double_array_descending(right_array, size - mid);

  // Merge the sorted halves.
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  while (i < mid && j < size - mid) {
    if (left_array[i] > right_array[j]) {
      array[k++] = left_array[i++];
    } else {
      array[k++] = right_array[j++];
    }
  }

  // Copy the remaining elements from the left half.
  while (i < mid) {
    array[k++] = left_array[i++];
  }

  // Copy the remaining elements from the right half.
  while (j < size - mid) {
    array[k++] = right_array[j++];
  }

  // Free the memory allocated for the left and right halves.
  free(left_array);
  free(right_array);
}
double randfrom(double min, double max)
{
  /*
  double mantissa = rand() % (1ull << 52) + 1;
  int exponent = rand() % (1024) - 1023;
  return mantissa * pow(2, exponent);
  */
  srand((unsigned)time(NULL));
  double X=((double)rand()/(double)RAND_MAX);
  return X;
}

void testtwosum(void (*implementation)(double, double, double *, double *) = twoSum)
{
  // test case modified as fast two sum gives sometimes results
  for (double b = 1; b < 100; b += 0.01)
  {
    double a = randfrom(-onedifference * onedifference, +onedifference * onedifference);
    double c = randfrom(-onedifference * onedifference, +onedifference * onedifference);

    double err_ours, err_ref;
    double res;

    implementation(a, c, &res, &err_ours);
    two_sum(a, c, &err_ref);
    assert(err_ref == err_ours || c == 0);
  }
}

void testfastfma(void (*implementation)(double, double, double *, double *) = twoMultFMA)
{
  // 2MultFMA
  for (double c = 1; c < 100; c += 0.01)
  {
    double a = randfrom(-onedifference, +onedifference);
    double b = randfrom(-onedifference, +onedifference);
    double res, err, error_ref;
    implementation(a, b, &res, &err);
    two_prod(a, b, &error_ref);
    assert(error_ref == err || b == 0);
  }
}

void testrenormalization(void (*implementation)(double *, int, double *, int) = renormalizationalgorithm)
{
  // Test size 1
  for (int c = 2; c < 100; c += 5)
  {
    double *renorm = new double[c];
    double *solution = new double[5];

    for (int i = 0; i < c; i++)
    {
      renorm[i] = allonesindouble;
    }
    implementation(renorm, c, solution, 1);
    renorm_rand2L(c, 1, renorm);
    for (int n = 0; n < 1; n++)
    {
      double rt = renorm[n];
      double st = solution[n];

      // test for reduction onto size 1 thus if failed thats the problem
      assert(abs(rt - st) < 0.00001);
    }

    delete[] renorm;
    delete[] solution;
  }

  for (int c = 2; c < 20; c += 1)
  {
    double *renorm = new double[c];
    double *solution = new double[5];

    for (int i = 0; i < c; i++)
    {
      renorm[i] = (allonesindouble * 1024) / (8 * i + 1);
    }
    implementation(renorm, c, solution, 2);
    renorm_rand2L(c, 2, renorm);

    for (int n = 0; n < 2; n++)
    {
      double rt = renorm[n];
      double st = solution[n];

      // test for reduction onto size 2 thus if failed thats the problem
      assert(abs(rt - st) < 0.00001);
    }

    delete[] renorm;
    delete[] solution;
  }

  for (int c = 3; c < 200; c += 1)
  {
    double *renorm = new double[c];
    double *solution = new double[5];

    for (int i = 0; i < c; i++)
    {
      renorm[i] = (allonesindouble * 1024) / (8 * i + 1);
    }
    implementation(renorm, c, solution, 3);

    renorm_rand2L(c, 3, renorm);
    for (int n = 0; n < 3; n++)
    {
      double rt = renorm[n];
      double st = solution[n];

      // test for reduction onto size 2 thus if failed thats the problem
      assert(abs(rt - st) < 0.00001);
    }

    delete[] renorm;
    delete[] solution;
  }

  // for (int c = 14; c < 50; c += 10)
  // {
  //   double *renorm = generate_d_nonoverlapping_expansion(c);
  //   double *solution = new double[c];

  //   implementation(renorm, c, solution, c);

  //   renorm_rand2L(c, c, renorm);
  //   for (int n = 0; n < c; n++)
  //   {
  //     double rt = renorm[n];
  //     double st = solution[n];

  //     // test for reduction onto size 2 thus if failed thats the problem
  //     assert(abs(rt - st) < 0.00001);
  //   }

  //   delete[] renorm;
  //   delete[] solution;
  // }
}

void testaddition(void (*implementation)(double *, double *, double *, int, int, int) = addition)
{
  // assumption that the numbers are still non p-2 overlapping simmilar to the renormalization
  // as addition is calling it and therefore it is kind of a natural implicit constraint
  for (int c = 2; c < 20; c += 1)
  {
    double *a = new double[c];
    double *b = new double[c];
    double *sol = new double[c]();
    double *sol_ref = new double[c]();
    for (int i = 0; i < c; i++)
    {
      a[i] = (allonesindouble * 1024) / (pow(8, 2 * i));
      b[i] = (allonesindouble * 1024) / (pow(8, 2 * i + 1));
    }
    // every 8 entries they go into the same "box"
    for (int i = 0; i < c; i += 8)
    {
      for (int x = 0; x < 8; x++)
      {
        if (i + x < c)
        {
          sol_ref[i / 8] += (a[i + x] + b[i + x]);
        }
      }
    }
    implementation(a, b, sol, c, c, c);

    for (int i = 0; i < c; i++)
    {
      double ref = sol_ref[i];
      double sol_our = sol[i];

      assert(abs(sol_our - ref) < 0.001);
    }

    delete[] a;
    delete[] b;
    delete[] sol;
    delete[] sol_ref;
  }
}

void testmultiplication(void (*implementation)(double *, double *, double *, int, int, int) = multiplication)
{

  // reference implementation certifiedMul
  // our implementation implementation
  // test case 1
  for (int c = 3; c < 10; c += 1)
  {
    double *a = new double[c];
    double *b = new double[c];
    double *sol = new double[c]();
    double *sol_ref = new double[c]();
    for (int i = 0; i < c; i++)
    {
      a[i] = (allonesindouble * 1024) / (pow(8, 2 * i));
      b[i] = (allonesindouble * 1024) / (pow(8, 2 * i + 1));
    }
    renormalizationalgorithm(a, c, a, c);
    renormalizationalgorithm(b, c, b, c);
    implementation(a, b, sol, c, c, 1);
    certifiedMul(c, c, 1, a, b, sol_ref);
    for (int i = 0; i < 1; i++)
    {
      double ref = sol_ref[i];
      double sol_our = sol[i];
      assert(abs(sol_our - ref) < 0.00001);
    }
    delete[] a;
    delete[] b;
    delete[] sol;
    delete[] sol_ref;
  }

  // test case 2
  for (int c = 3; c < 15; c += 1)
  {
    double *a = new double[c];
    double *b = new double[c];
    double *sol = new double[c]();
    double *sol_ref = new double[c]();
    for (int i = 0; i < c; i++)
    {
      a[i] = (allonesindouble * 1024) / (pow(8, 2 * i));
      b[i] = (allonesindouble * 1024) / (pow(8, 2 * i + 1));
    }

    renormalizationalgorithm(a, c, a, c);
    renormalizationalgorithm(b, c, b, c);
    implementation(a, b, sol, c, c, 2);
    certifiedMul(c, c, 2, a, b, sol_ref);
    for (int i = 0; i < 2; i++)
    {
      double ref = sol_ref[i];
      double sol_our = sol[i];
      assert(abs(sol_our - ref) < 0.00000001);
    }
    delete[] a;
    delete[] b;
    delete[] sol;
    delete[] sol_ref;
  }
  // test case 3
  for (int c = 3; c < 15; c += 1)
  {
    double *a = new double[c];
    double *b = new double[c];
    double *sol = new double[c]();
    double *sol_ref = new double[c]();
    for (int i = 0; i < c; i++)
    {
      a[i] = (allonesindouble * 1024) / (pow(8, 2 * i));
      b[i] = (allonesindouble * 1024) / (pow(8, 2 * i + 1));
    }
    renormalizationalgorithm(a, c, a, c);
    renormalizationalgorithm(b, c, b, c);
    implementation(a, b, sol, c, c, 3);
    certifiedMul(c, c, 3, a, b, sol_ref);
    for (int i = 0; i < 3; i++)
    {
      double ref = sol_ref[i];
      double sol_our = sol[i];
      assert(abs(sol_our - ref) < 0.00000001);
    }
    delete[] a;
    delete[] b;
    delete[] sol;
    delete[] sol_ref;
  }
  for(int p = 3; p < 400; p += 1){
    for (int c = 3; c < 15; c += 1)
    {
      double *a = new double[c];
      double *b = new double[c];
      double *sol = new double[c]();
      double *sol_ref = new double[c]();
      for (int i = 0; i < c; i++)
      {
        a[i] = randfrom(0,1);
        b[i] =  randfrom(0,1);
      }
      renormalizationalgorithm(a, c, a, c);
      renormalizationalgorithm(b, c, b, c);
      // as the mul alorithm implicitly assumes that the input is sorted sort it:
      sort_double_array_descending(a, c);
      sort_double_array_descending(b, c);
      implementation(a, b, sol, c, c, 3);
      certifiedMul(c, c, 3, a, b, sol_ref);
      for (int i = 0; i < 3; i++)
      {
        double ref = sol_ref[i];
        double sol_our = sol[i];
        assert(abs(sol_our - ref) <= (abs(0.000000001* a[0]*b[0])));
      }
      delete[] a;
      delete[] b;
      delete[] sol;
      delete[] sol_ref;
    }
  }
}

void testfourmultiplication(void (*implementation)(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, int, int) = fourtimesmultiplicationversion3)
{

  for (int c = 3; c < 5; c += 1)
  {
    double *a0 = new double[c];
    double *b0 = new double[c];
    double *sol0 = new double[4];
    double *sol_ref0 = new double[4];
    double *a1 = new double[c];
    double *b1 = new double[c];
    double *sol1 = new double[4];
    double *sol_ref1 = new double[4];
    double *a2 = new double[c];
    double *b2 = new double[c];
    double *sol2 = new double[4];
    double *sol_ref2 = new double[4];
    double *a3 = new double[c];
    double *b3 = new double[c];
    double *sol3 = new double[4];
    double *sol_ref3 = new double[4];

    for (int i = 0; i < c; i++)
    {
      a0[i] = randfrom(0,1);
      b0[i] = randfrom(0,1);
      a1[i] =randfrom(0,1);
      b1[i] = randfrom(0,1);
      a2[i] = randfrom(0,1);
      b2[i] = randfrom(0,1);
      a3[i] = randfrom(0,1);
      b3[i] = randfrom(0,1);
    }
    
    renormalizationalgorithm(a0, c, a0, c);
    renormalizationalgorithm(b0, c, b0, c);
    renormalizationalgorithm(a1, c, a1, c);
    renormalizationalgorithm(b1, c, b1, c);
    renormalizationalgorithm(a2, c, a2, c);
    renormalizationalgorithm(b2, c, b2, c);
    renormalizationalgorithm(a3, c, a3, c);
    renormalizationalgorithm(b3, c, b3, c);
    
    implementation(a0, b0, a1, b1, a2, b2, a3, b3, sol0, sol1, sol2, sol3, c, c, 4);
    
    multiplication(a0, b0, sol_ref0, c, c,4);
    multiplication(a1, b1, sol_ref1, c, c,4);
    multiplication(a2, b2, sol_ref2, c, c,4);
    multiplication(a3, b3, sol_ref3, c, c,4);
 
    
    for (int i = 0; i < 4; i++)
    {
      double ref0 = sol_ref0[i];
      double sol_our0 = sol0[i];
      double ref1 = sol_ref1[i];
      double sol_our1 = sol1[i];
      double ref2 = sol_ref2[i];
      double sol_our2 = sol2[i];
      double ref3 = sol_ref3[i];
      double sol_our3 = sol3[i];
      assert(abs(sol_our0 - ref0) == (0));
      assert(abs(sol_our1 - ref1) == (0));
      assert(abs(sol_our2 - ref2) == (0));
      assert(abs(sol_our3 - ref3) == (0));
     
    }

    delete[] a0;
    delete[] b0;
    delete[] sol0;
    delete[] sol_ref0;
    delete[] a1;
    delete[] b1;
    delete[] sol1;
    delete[] sol_ref1;
    delete[] a2;
    delete[] b2;
    delete[] sol2;
    delete[] sol_ref2;
    delete[] a3;
    delete[] b3;
    delete[] sol3;
    delete[] sol_ref3;
  }
}