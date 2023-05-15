#define NOMINMAX

#include <limits>
#include <random>

#include "basefunctions.h"

#include "mult2_optimizations.cpp"
#include "timing.cpp"

// compile with g++ -std=c++17 ./benchmark.cpp
// TODO: make sure input is non-S-overlapping, check validity of d-non-overlapping expansion
//       calculate median of measurements, remove comment about validity of non-overlap
//       take out memory allocation from measured functions

#define DEFAULT_TERM_COUNT 256

double generate_random_double()
{
    random_device random_device;
    mt19937 generator(random_device());

    double double_minimum = numeric_limits<double>::min();
    double double_maximum = numeric_limits<double>::max();

    default_random_engine random_engine;
    uniform_real_distribution<double> distribution(0., 1.);

    return distribution(generator);
}

// `term_count` needs to be at most 1022
// Note: this function does not necessarily generate a valid d-digit-non-overlapping representation
// (inequality (4) in the definition is not always met), it's only used for benchmarking
double *generate_maximally_overlapping_expansion(int term_count)
{
    int exponent = 1023;
    double mantissa = generate_random_double();

    double *terms = new double[term_count];

    for (int i = 0; i < term_count; i++)
    {
        terms[i] = ldexp(mantissa, exponent - 2 * i);
    }

    return terms;
}

double *generate_ulp_nonoverlapping_expansion(int term_count)
{
    int exponent = 500;
    double mantissa = generate_random_double();

    double *terms = new double[term_count];

    for (int i = 0; i < term_count; i++)
    {
        terms[i] = ldexp(mantissa, exponent);
        exponent -= 53;
    }

    return terms;
}

template <typename T_compute_return, class T_compute_arguments>
vector<double> run_measurements(
    double (*measure)(T_compute_return (*implementation)(T_compute_arguments...), int term_count),
    T_compute_return (*implementation)(T_compute_arguments...))
{
    vector<double> runtimes = {};

    for (int term_count = 22; term_count < 1023; term_count += 50)
    {
        runtimes.push_back(measure(implementation, term_count));
    }

    return runtimes;
}

// Measurement functions

double measure_two_sum(void (*implementation)(double, double, double *, double *) = twoSum)
{
    double a = generate_random_double();
    double b = generate_random_double();

    double *result = new double;
    double *error = new double;

    double runtime = measure_runtime(implementation, a, b, result, error);

    delete result;
    delete error;

    return runtime;
}

double measure_two_mult_FMA(void (*implementation)(double, double, double *, double *) = twoMultFMA)
{
    double a = generate_random_double();
    double b = generate_random_double();

    double *result = new double;
    double *error = new double;

    double runtime = measure_runtime(implementation, a, b, result, error);

    delete result;
    delete error;

    return runtime;
}

double measure_vec_sum(
    void (*implementation)(double *, double *, int) = vecSum,
    int term_count = DEFAULT_TERM_COUNT)
{
    double *expansion = generate_maximally_overlapping_expansion(term_count);
    double *result = new double[term_count];

    double runtime = measure_runtime(implementation, expansion, result, term_count);

    delete[] expansion;
    delete[] result;

    return runtime;
}

double measure_vec_sum_err_branch(
    void (*implementation)(double *, int, int, double *) = vecSumErrBranch,
    int term_count = DEFAULT_TERM_COUNT)
{
    double *expansion = generate_maximally_overlapping_expansion(term_count);
    double *result = new double[term_count];

    double runtime = measure_runtime(implementation, expansion, term_count, term_count,
                                     result);

    delete[] expansion;
    delete[] result;

    return runtime;
}

double measure_vec_sum_err(
    void (*implementation)(double *, int, double *) = vecSumErr,
    int term_count = DEFAULT_TERM_COUNT)
{
    double *expansion = generate_maximally_overlapping_expansion(term_count);
    double *result = new double[term_count];

    double runtime = measure_runtime(implementation, expansion, term_count, result);

    delete[] expansion;
    delete[] result;

    return runtime;
}

double measure_renormalization(
    void (*implementation)(double *, int, double *, int) = renormalizationalgorithm,
    int term_count = DEFAULT_TERM_COUNT)
{
    double *expansion = generate_maximally_overlapping_expansion(term_count);
    double *result = generate_maximally_overlapping_expansion(term_count);

    double runtime = measure_runtime(implementation, expansion, term_count, result, term_count);

    delete[] expansion;
    delete[] result;

    return runtime;
}

double measure_addition(
    void (*implementation)(double *, double *, double *, int, int, int) = addition,
    int term_count = DEFAULT_TERM_COUNT)
{
    double *a = generate_maximally_overlapping_expansion(term_count);
    double *b = generate_maximally_overlapping_expansion(term_count);

    double *result = new double[term_count];

    double runtime = measure_runtime(implementation, a, b, result, term_count, term_count,
                                     term_count);

    delete[] a;
    delete[] b;
    delete[] result;

    return runtime;
}

double measure_multiplication(
    void (*implementation)(double *, double *, double *, int, int, int) = multiplication,
    int term_count = DEFAULT_TERM_COUNT)
{
    double *a = generate_maximally_overlapping_expansion(term_count);
    double *b = generate_maximally_overlapping_expansion(term_count);

    double *result = new double[term_count];

    double runtime = measure_runtime(implementation, a, b, result, term_count, term_count,
                                     term_count);

    delete[] a;
    delete[] b;
    delete[] result;

    return runtime;
}

double measure_multiplication2(
    void (*implementation)(double *, double *, double *, int, int, int) = mult2,
    int term_count = DEFAULT_TERM_COUNT)
{
    double *a = generate_ulp_nonoverlapping_expansion(term_count);
    double *b = generate_ulp_nonoverlapping_expansion(term_count);

    double *result = new double[term_count];

    double runtime = measure_runtime(implementation, a, b, result, term_count, term_count, term_count);

    delete[] a;
    delete[] b;
    delete[] result;

    return runtime;
}
