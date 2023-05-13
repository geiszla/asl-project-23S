#include <filesystem>
#include <fstream>
#include <iostream>

#include "basefunctions.cpp"
#include "timing.cpp"

// compile with g++ -std=c++17 ./benchmark.cpp
// TODO: make sure input is non-overlapping, calculate median of measurements,
//       clean up `pyproject.toml` (authors, repo link, license. etc.)

#define OUTPUT_PATH "../results"

// Flop counts

#define TWO_SUM_FLOPS 6
#define TWO_MULT_FMA_FLOPS 2

unsigned int get_vec_sum_flops(int vector_length)
{
    return (vector_length - 1) * TWO_SUM_FLOPS;
}

unsigned int get_vec_sum_err_branch_flops(int input_expansion_length)
{
    return (input_expansion_length - 1) * TWO_SUM_FLOPS;
}

unsigned int get_vec_sum_err_flops(int expansion_length)
{
    return (expansion_length - 1) * TWO_SUM_FLOPS;
}

unsigned int get_renormalization_flops(int input_expansion_length, int output_expansion_length)
{
    int vec_sum_flops = get_vec_sum_flops(input_expansion_length);
    int vec_sum_err_branch_flops = get_vec_sum_err_branch_flops(input_expansion_length);
    int vec_sum_err_flops = get_vec_sum_err_flops(
        (output_expansion_length * (output_expansion_length + 1) / 2.0));

    return vec_sum_flops + vec_sum_err_branch_flops + vec_sum_err_flops;
}

unsigned int get_addition_flops(int a_length, int b_length, int result_length)
{
    return get_renormalization_flops(a_length + b_length, result_length);
}

unsigned int get_multiplication_flops(int expansion_length)
{
    return expansion_length * (expansion_length + 1) / 2 * TWO_MULT_FMA_FLOPS +
           (expansion_length / 2) * get_vec_sum_flops(expansion_length) +
           (expansion_length - 1) * 2 +
           pow(expansion_length, 2) +
           get_renormalization_flops(expansion_length + 1, expansion_length);
}

// Helpers

double generate_random_double(double max)
{
    srand(time(nullptr));
    return (max / RAND_MAX) * rand();
}

void fill_vector(double x[], size_t size)
{
    // -1022 to +1023 in double
    int exponent = 500;

    for (size_t i = 0; i < size; i++)
    {
        x[i] = ldexp(generate_random_double(1), exponent);
    }
}

void write_results(string algorithm_name, double runtime, double performance, int input_size,
                   int input_size_exponent, ostream &output_file)
{
    output_file << algorithm_name << ";" << input_size << ";" << runtime << ";" << performance << endl;

    cout << "Input size: 2^" << input_size_exponent
         << "     runtime: " << runtime << " cycles, performance: "
         << performance << " flops/cycles" << endl;
}

// Benchmarking functions

void benchmark_two_sum()
{
    double a = generate_random_double(1);
    double b = generate_random_double(1);

    double *result = new double;
    double *error = new double;

    double runtime = measure_runtime(&twoSum, a, b, result, error);

    delete result;
    delete error;

    cout << "2Sum\t\truntime: " << runtime << " cycles, performance: "
         << (TWO_SUM_FLOPS / runtime) << " flops/cycles" << endl;
}

void benchmark_two_mult_FMA()
{
    double a = generate_random_double(1);
    double b = generate_random_double(1);

    double *result = new double;
    double *error = new double;

    double runtime = measure_runtime(&twoMultFMA, a, b, result, error);

    delete result;
    delete error;

    cout << "2MultFMA\truntime: " << runtime << " cycles, performance: "
         << (TWO_MULT_FMA_FLOPS / runtime) << " flops/cycles" << endl;
}

void benchmark_vec_sum(ofstream &output_file)
{
    cout << endl
         << "VecSum: " << endl;

    for (size_t i = 4; i < 15; i++)
    {
        int input_size = pow(2, i);

        double *x = new double[input_size];
        fill_vector(x, input_size);

        double *errors = new double[input_size];

        double runtime = measure_runtime(&vecSum, x, errors, input_size);
        double performance = get_vec_sum_flops(input_size) / runtime;

        write_results("VecSum", runtime, performance, input_size, i, output_file);
    }
}

void benchmark_vec_sum_err_branch(ofstream &output_file)
{
    cout << endl
         << "VecSumErrBranch: " << endl;

    for (size_t i = 4; i < 15; i++)
    {
        int input_size = pow(2, i);

        double *expansion = new double[input_size];
        fill_vector(expansion, input_size);

        double runtime = measure_runtime(&vecSumErrBranch, expansion, input_size, input_size);
        double performance = get_vec_sum_err_branch_flops(input_size) / runtime;

        write_results("VecSumErrBranch", runtime, performance, input_size, i, output_file);
    }
}

void benchmark_vec_sum_err(ofstream &output_file)
{
    cout << endl
         << "VecSumErr: " << endl;

    for (size_t i = 4; i < 15; i++)
    {
        int input_size = pow(2, i);

        double *expansion = new double[input_size];
        fill_vector(expansion, input_size);

        double runtime = measure_runtime(&vecSumErr, expansion, input_size);
        double performance = get_vec_sum_err_flops(input_size) / runtime;

        write_results("VecSumErr", runtime, performance, input_size, i, output_file);
    }
}

void benchmark_renormalization(ofstream &output_file)
{
    cout << endl
         << "Renormalization: " << endl;

    for (size_t i = 4; i < 15; i++)
    {
        int input_size = pow(2, i);

        double *x = new double[input_size];
        fill_vector(x, input_size);

        double *expansion = new double[input_size];
        fill_vector(expansion, input_size);

        double runtime = measure_runtime(
            &renormalizationalgorithm, x, input_size, expansion, input_size);
        double performance = get_renormalization_flops(input_size, input_size) / runtime;

        write_results("Renormalization", runtime, performance, input_size, i, output_file);
    }
}

void benchmark_addition(ofstream &output_file)
{
    cout << endl
         << "Addition: " << endl;

    for (size_t i = 4; i < 20; i++)
    {
        int input_size = pow(2, i);

        double *a = new double[input_size];
        fill_vector(a, input_size);

        double *b = new double[input_size];
        fill_vector(b, input_size);

        double *result = new double[input_size];

        double runtime = measure_runtime(&addition, a, b, result, input_size, input_size,
                                         input_size);

        double performance = get_addition_flops(input_size, input_size, input_size) / runtime;

        write_results("Addition", runtime, performance, input_size, i, output_file);
    }
}

void benchmark_multiplication(ofstream &output_file)
{
    cout << endl
         << "Multiplication: " << endl;

    for (size_t i = 4; i < 20; i++)
    {
        int input_size = pow(2, i);

        double *a = new double[input_size];
        fill_vector(a, input_size);

        double *b = new double[input_size];
        fill_vector(b, input_size);

        double *result = new double[input_size];

        double runtime = measure_runtime(&multiplication, a, b, result);
        double performance = get_multiplication_flops(input_size) / runtime;

        write_results("Multiplication", runtime, performance, input_size, i, output_file);
    }
}

// Main

int main()
{
    error_code directory_error;

    // Create directory if does not exist
    if ((!filesystem::is_directory(OUTPUT_PATH) || !filesystem::exists(OUTPUT_PATH)) &&
        !filesystem::create_directory(OUTPUT_PATH, directory_error))
    {
        cerr << "Error: cannot create output directory: " << directory_error.message() << endl;
        return 1;
    }

    // Open output file
    ofstream output_file;
    output_file.open(string(OUTPUT_PATH) + "/benchmark.csv");

    if (!output_file.is_open())
    {
        cerr << "Error: cannot open output file: " << strerror(errno) << endl;
        return 1;
    }

    output_file << "Algorithm;Input size;Runtime;Performance" << endl;

    cout.setf(ios::fixed);
    cout.precision(2);

    benchmark_two_sum();
    benchmark_two_mult_FMA();

    benchmark_vec_sum(output_file);
    benchmark_vec_sum_err_branch(output_file);
    benchmark_vec_sum_err(output_file);

    benchmark_renormalization(output_file);
    benchmark_addition(output_file);
    benchmark_multiplication(output_file);

    output_file.close();
}
