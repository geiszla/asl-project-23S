#define NOMINMAX

#include <limits>
#include <random>

#include "timing.cpp"

// Needs to be low, otherwise we cannot place expansions on the stack
// Needs to be at most 29 to comply with ulp-nonoverlapping restrictions
#define DEFAULT_TERM_COUNT 29

double generate_random_double()
{
    random_device random_device;
    mt19937 generator(random_device());

    double double_minimum = numeric_limits<double>::min();
    double double_maximum = numeric_limits<double>::max();

    default_random_engine random_engine;
    uniform_real_distribution<double> double_distribution(1., 2.);

    return double_distribution(generator);
}

double generate_random_mantissa(int ending_zeros_count)
{
    random_device random_device;
    mt19937 generator(random_device());
    uniform_int_distribution<mt19937::result_type> boolean_distribution(0, 1);

    double mantissa = 1;

    for (int i = 1; i < 53 - ending_zeros_count; i++)
    {
        mantissa += boolean_distribution(generator) * pow(2, -i);
    }

    return mantissa;
}

/**
 * Generates an expansion that is at most d-overlapping (with d = 2)
 * @param term_count The number of terms to generate. Needs to be at most 742
 */
double *generate_d_nonoverlapping_expansion(int term_count)
{
    int exponent = 511;

    double *terms = new double[term_count];

    terms[0] = ldexp(generate_random_double(), exponent);

    for (int i = 1; i < term_count; i++)
    {
        // For this to be at most d-overlapping, exponent needs to decrease by at least 2
        // term-by-term and there needs to be at least 49 zeros at the end of each term
        terms[i] = ldexp(generate_random_mantissa(49), exponent - 51 - 2 * (i - 1));
    }

    return terms;
}

/**
 * Generates an expansion that is S-nonoverlapping
 * @param term_count The number of terms to generate. Needs to be at most 1533
 */
double *generate_s_nonoverlapping_expansion(int term_count)
{
    // Generate `term_count` number of terms wit the minimum number of ending zeros possible
    int exponent_difference = 1533 / (term_count - 1);
    int ending_zeros_count = max(0, 53 - exponent_difference);

    int exponent = 511;

    double *terms = new double[term_count];

    for (int i = 0; i < term_count; i++)
    {
        // For this to be at most d-overlapping, exponent needs to decrease by at least 2
        // term-by-term and there needs to be at least 49 zeros at the end of each term
        terms[i] = ldexp(generate_random_mantissa(ending_zeros_count), exponent);
        exponent -= exponent_difference;
    }

    return terms;
}

/**
 * Generates an expansion that is ulp-overlapping
 * @param term_count The number of terms to generate. Needs to be at most 29
 */
double *generate_ulp_nonoverlapping_expansion(int term_count)
{
    int exponent = 511;
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
    double *expansion = generate_d_nonoverlapping_expansion(term_count);
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
    double *expansion = generate_s_nonoverlapping_expansion(term_count);
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
    double *expansion = generate_d_nonoverlapping_expansion(term_count);
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
    double *expansion = generate_d_nonoverlapping_expansion(term_count);
    double *result = generate_d_nonoverlapping_expansion(term_count);

    double runtime = measure_runtime(implementation, expansion, term_count, result, term_count);

    delete[] expansion;
    delete[] result;

    return runtime;
}

double measure_addition(
    void (*implementation)(double *, double *, double *, int, int, int) = addition,
    int term_count = DEFAULT_TERM_COUNT)
{
    double *a = generate_d_nonoverlapping_expansion(term_count);
    double *b = generate_d_nonoverlapping_expansion(term_count);

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
    double *a = generate_d_nonoverlapping_expansion(term_count);
    double *b = generate_d_nonoverlapping_expansion(term_count);

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
