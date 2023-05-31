#include <iostream>
#include <map>

#include "measure.cpp"
#include "testalgo5and6.cpp"
#include "testAlgos7and8.cpp"

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

template <typename T_test_function>
map<Algorithm, T_test_function> tests{
    {two_sum_algorithm, (T_test_function)testtwosum},
    {vec_sum_err_algorithm, (T_test_function)testvecsumerr},
    {vec_sum_err_branch_algorithm, (T_test_function)testvecsumerrbranch},
    {two_mult_FMA_algorithm, (T_test_function)testfastfma},
    {renormalization_algorithm, (T_test_function)testrenormalization},
    {addition_algorithm, (T_test_function)testaddition},
    {multiplication_algorithm, (T_test_function)testmultiplication}};

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
 * @param implementation The implementated function to evaluate
 * @param name Name of the implementation (only used to display in command line)
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
        cout << "; " <<  naive_runtime / runtime << "x speedup" <<endl;
    }
    else
    {
        cout << endl;
    }
}
