#include "evaluate.cpp"

extern "C"
{
  #include "vec_sum_optimizations.c"
}

// Compile with `g++ -std=c++17 ./benchmark.cpp` (plus additional optimization flags)

int main()
{
  // Pass the naive implementation's runtime to `evaluate_implementation` to calculate speedup
  double naive_runtime = measure_vec_sum();

  evaluate_implementation(renormalizationalgorithm2, "renormalization2", renormalization_algorithm);
  evaluate_implementation(vecSum3, "vecSum3", vec_sum_algorithm, naive_runtime);
}
