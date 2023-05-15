#include <iostream>
#include <map>

#include "measure.cpp"
#include "testalgo5and6.cpp"

enum Algorithm
{
    two_sum_algorithm,
    two_mult_FMA_algorithm,
    vec_sum_algorithm,
    vec_sum_err_branch_algorithm,
    vec_sum_err_algorithm,
    renormalization_algorithm,
    addition_algorithm,
    multiplication_algorithm,
    multiplication2_algorithm
};

template <typename T_compute_function>
map<Algorithm, T_compute_function> tests{
    {two_sum_algorithm, (T_compute_function)testtwosum},
    {two_mult_FMA_algorithm, (T_compute_function)testfastfma},
    {renormalization_algorithm, (T_compute_function)testrenormalization},
    {addition_algorithm, (T_compute_function)testaddition}};

template <typename T_measure_function>
map<Algorithm, T_measure_function> measurement_functions{
    {two_sum_algorithm, (T_measure_function)measure_two_sum},
    {two_mult_FMA_algorithm, (T_measure_function)measure_two_mult_FMA},
    {vec_sum_algorithm, (T_measure_function)measure_vec_sum},
    {vec_sum_err_branch_algorithm, (T_measure_function)measure_vec_sum_err_branch},
    {vec_sum_err_algorithm, (T_measure_function)measure_vec_sum_err},
    {renormalization_algorithm, (T_measure_function)measure_renormalization},
    {addition_algorithm, (T_measure_function)measure_addition},
    {multiplication_algorithm, (T_measure_function)measure_multiplication},
    {multiplication2_algorithm, (T_measure_function)measure_multiplication2}};

/**
 * Evaluates an implementation of an algorithm
 *
 * @param implementation The function implementation to evaluate
 * @param name Name of the implementation (only used for display in command line)
 * @param algorithm The algorithm from `Algorithm` enum that is being implemented
 * @param naive_runtime If passed, the speedup of the implementation compared to this
 * will be calculated
 */
template <typename T_compute_return, class... T_compute_arguments>
void evaluate_implementation(
    T_compute_return (*implementation)(T_compute_arguments...),
    string name,
    Algorithm algorithm,
    double naive_runtime = 0.)
{
    auto typed_tests = tests<void (*)(T_compute_return (*)(T_compute_arguments...))>;

    if (typed_tests.count(algorithm) > 0)
    {
        typed_tests[algorithm](implementation);
    }
    else
    {
        cerr << "No test were found for \"" << name << '"' << endl;
    }

    auto measurement_function = measurement_functions<
        double (*)(T_compute_return (*)(T_compute_arguments...), int)>[algorithm];

    double runtime = measurement_function(implementation, DEFAULT_TERM_COUNT);

    cout.setf(ios::fixed);
    cout.precision(3);

    cout << name << " - Runtime: " << runtime;

    if (naive_runtime != 0)
    {
        cout << "; " << runtime / naive_runtime << "x speedup";
    }
    else
    {
        cout << endl;
    }
}
