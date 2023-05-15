#include <iostream>
#include <map>

#include "measure.cpp"
#include "testalgo5and6.cpp"

enum FPFunction
{
    two_sum_function,
    two_mult_FMA_function,
    vec_sum_function,
    vec_sum_err_branch_function,
    vec_sum_err_function,
    renormalization_function,
    addition_function,
    multiplication_function,
    multiplication2_function
};

template <typename T_compute_function>
map<FPFunction, T_compute_function> tests{
    {two_sum_function, (T_compute_function)testtwosum},
    {two_mult_FMA_function, (T_compute_function)testfastfma},
    {renormalization_function, (T_compute_function)testrenormalization},
    {addition_function, (T_compute_function)testaddition}};

template <typename T_measure_function>
map<FPFunction, T_measure_function> measurement_functions{
    {two_sum_function, (T_measure_function)measure_two_sum},
    {two_mult_FMA_function, (T_measure_function)measure_two_mult_FMA},
    {vec_sum_function, (T_measure_function)measure_vec_sum},
    {vec_sum_err_branch_function, (T_measure_function)measure_vec_sum_err_branch},
    {vec_sum_err_function, (T_measure_function)measure_vec_sum_err},
    {renormalization_function, (T_measure_function)measure_renormalization},
    {addition_function, (T_measure_function)measure_addition},
    {multiplication_function, (T_measure_function)measure_multiplication},
    {multiplication2_function, (T_measure_function)measure_multiplication2}};

template <typename T_compute_return, class... T_compute_arguments>
void evaluate_implementation(
    T_compute_return (*implementation)(T_compute_arguments...),
    string name,
    FPFunction function_type,
    double naive_runtime = 0.)
{
    auto typed_tests = tests<void (*)(T_compute_return (*)(T_compute_arguments...))>;

    if (typed_tests.count(function_type) > 0)
    {
        typed_tests[function_type](implementation);
    }
    else
    {
        cerr << "No test were found for \"" << name << '"' << endl;
    }

    auto measurement_function = measurement_functions<
        double (*)(T_compute_return (*)(T_compute_arguments...), int)>[function_type];

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
