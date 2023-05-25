#include "evaluate.cpp"

extern "C"
{
  #include "vec_sum_optimizations.c"
  #include "vec_sum_err_optimizations.c"
  #include "vec_sum_err_branch_optimizations.c"
  #include "renormalization_optimizations.c"
  #include "addition_optimizations.c"
}

// Compile with `g++ -std=c++17 ./benchmark.cpp` (plus additional optimization flags)

int main()
{
  /* vecSum */
  // Pass the naive implementation's runtime to `evaluate_implementation` to calculate speedup
  // double naive_runtime = measure_vec_sum();

  // evaluate_implementation(renormalizationalgorithm2, "renormalization2", renormalization_algorithm);
  // evaluate_implementation(vecSum3, "vecSum3", vec_sum_algorithm, naive_runtime);

  /* vecSumErr */
  // double naive_runtime = measure_vec_sum_err();

  // evaluate_implementation(vecSumErr2, "vecSumErr2", vec_sum_err_algorithm, naive_runtime);

  /* vecSumErrBranch */
  // double naive_runtime = measure_vec_sum_err_branch();

  // evaluate_implementation(vecSumErrBranch2, "vecSumErrBranch2", vec_sum_err_branch_algorithm,
  //                         naive_runtime);

  /* renormalization */
  // double naive_runtime = measure_renormalization();

  // evaluate_implementation(renormalization3, "renormalization3", renormalization_algorithm,
  //                         naive_runtime);
  // evaluate_implementation(renormalization4, "renormalization4", renormalization_algorithm,
  //                         naive_runtime);

  /* addition */
  double naive_runtime = measure_addition();

  evaluate_implementation(addition3, "addition3", addition_algorithm, naive_runtime);
  evaluate_implementation(addition4, "addition4", addition_algorithm, naive_runtime);
}
