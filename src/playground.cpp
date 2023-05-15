#include "evaluate.cpp"

// Compile with `g++ -std=c++17 ./benchmark.cpp` (plus additional optimization flags)

int main()
{
  // Pass the naive implementation's runtime to `evaluate_implementation` to calculate speedup
  double naive_runtime = measure_vec_sum();

  evaluate_implementation(vecSum, "vec_sum_slow", vec_sum_algorithm, naive_runtime);
}
