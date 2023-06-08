extern "C"
{
#include "reference.h"
#include "basefunctions.c"
#include "vec_sum_optimizations.c"
#include "vec_sum_err_optimizations.c"
#include "vec_sum_err_branch_optimizations.c"
#include "renormalization_optimizations.c"
#include "addition_optimizations.c"
#include "multiplication_optimizations.c"
}

#include "evaluate.cpp"

// Compile with `g++ -std=c++17 ./benchmark.cpp` (plus additional optimization flags)

int main()
{
  /* vecSum */
  // Pass the naive implementation's runtime to `evaluate_implementation` to calculate speedup
  // double naive_runtime = measure_vec_sum();

  // evaluate_implementation(renormalizationalgorithm2, "renormalization2",
  //                         renormalization_algorithm);
  // evaluate_implementation(vecSum6, "vecSum6", vec_sum_algorithm, naive_runtime);
  // evaluate_implementation(vecSum6_fast, "vecSum6_fast", vec_sum_algorithm, naive_runtime);

  /* vecSumErrBranch */
  // double naive_runtime = measure_vec_sum_err_branch();

  // evaluate_implementation(vecSumErrBranch2, "vecSumErrBranch2", vec_sum_err_branch_algorithm,
  //                         naive_runtime);
  // evaluate_implementation(vecSumErrBranch_fast, "vecSumErrBranch_fast",
  //                         vec_sum_err_branch_algorithm, naive_runtime);

  /* vecSumErr */
  // double naive_runtime = measure_vec_sum_err();

  // evaluate_implementation(vecSumErr2, "vecSumErr2", vec_sum_err_algorithm, naive_runtime);
  // evaluate_implementation(vecSumErr_fast, "vecSumErr_fast", vec_sum_err_algorithm, naive_runtime);

  /* renormalization */
  // double naive_runtime = measure_renormalization();

  // evaluate_implementation(renormalization4, "renormalization4", renormalization_algorithm,
  //                         naive_runtime);
  // evaluate_implementation(renormalization_fast, "renormalization_fast", renormalization_algorithm,
  //                         naive_runtime);

  /* addition */
  // double naive_runtime = measure_addition();

  // evaluate_implementation(addition4, "addition4", addition_algorithm, naive_runtime);
  // evaluate_implementation(addition_fast, "addition_fast", addition_algorithm, naive_runtime);

  /* multiplication */
  double naive_runtime = measure_multiplication();

  evaluate_implementation(multiplication2, "multiplication2", multiplication_algorithm,
                          naive_runtime);
  evaluate_implementation(multiplication3, "multiplication3", multiplication_algorithm,
                          naive_runtime);
}
