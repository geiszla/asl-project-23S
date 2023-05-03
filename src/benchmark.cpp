#include <filesystem>
#include <fstream>
#include <iostream>

#include "basefunctions.cpp"
#include "timing.cpp"

// compile with g++ -std=c++17 ./benchmark.cpp
// TODO: measure all functions only once, factor out common code

#define OUTPUT_PATH "../results"

// Helpers

double generate_random_double(double max)
{
    srand(time(nullptr));
    return (max / RAND_MAX) * rand();
}

void fill_vector(double x[], size_t size)
{
    // -1022 to +1023 in double
    int exponent = 1023;

    for (size_t i = 0; i < size; i++)
    {
        x[i] = ldexp(generate_random_double(1), exponent);
    }
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

    cout << "2Sum runtime: " << runtime << " cycles" << endl;
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

    cout << "2MultFMA runtime: " << runtime << " cycles" << endl;
}

void benchmark_vec_sum(ofstream &output_file)
{
    cout << endl << "VecSum runtime: " << endl;

    for (size_t i = 4; i < 20; i++)
    {
        int input_size = pow(2, i);

        double *x = new double[input_size];
        fill_vector(x, input_size);

        double *errors = new double[input_size];

        double runtime = measure_runtime(&vecSum, x, errors);

        output_file << input_size << ";" << runtime << endl;
        cout << "Input size: 2^" << i << ", Runtime: " << runtime << " cycles" << endl;
    }
}

void benchmark_vec_sum_err_branch(ofstream &output_file)
{
    cout << endl << "VecSumErrBranch runtime: " << endl;

    for (size_t i = 4; i < 20; i++)
    {
        int input_size = pow(2, i);

        double *expansion = new double[input_size];
        fill_vector(expansion, input_size);

        double runtime = measure_runtime(&vecSumErrBranch, expansion, input_size, input_size);

        output_file << input_size << ";" << runtime << endl;
        cout << "Input size: 2^" << i << ", Runtime: " << runtime << " cycles" << endl;
    }
}

void benchmark_vec_sum_err(ofstream &output_file)
{
    cout << endl << "VecSumErr runtime: " << endl;

    for (size_t i = 4; i < 20; i++)
    {
        int input_size = pow(2, i);

        double *expansion = new double[input_size];
        fill_vector(expansion, input_size);

        double runtime = measure_runtime(&vecSumErr, expansion, input_size);

        output_file << input_size << ";" << runtime << endl;
        cout << "Input size: 2^" << i << ", Runtime: " << runtime << " cycles" << endl;
    }
}

void benchmark_renormalization(ofstream &output_file)
{
    cout << endl << "Renormalization runtime: " << endl;

    for (size_t i = 4; i < 20; i++)
    {
        int input_size = pow(2, i);

        double *x = new double[input_size];
        fill_vector(x, input_size);

        double *expansion = new double[input_size];
        fill_vector(expansion, input_size);

        double runtime = measure_runtime(&renormalizationalgorithm, x, expansion, input_size);

        output_file << input_size << ";" << runtime << endl;
        cout << "Input size: 2^" << i << ", Runtime: " << runtime << " cycles" << endl;
    }
}

void benchmark_multiplication(ofstream &output_file)
{
    cout << endl << "Multiplication runtime: " << endl;

    for (size_t i = 4; i < 20; i++)
    {
        int input_size = pow(2, i);

        double *a = new double[input_size];
        fill_vector(a, input_size);

        double *b = new double[input_size];
        fill_vector(b, input_size);

        double *result = new double[input_size];

        double runtime = measure_runtime(&multiplication, a, b, result);

        output_file << input_size << ";" << runtime << endl;
        cout << "Input size: 2^" << i << ", Runtime: " << runtime << " cycles" << endl;
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

    output_file << "Input size;Runtime" << endl;

    cout.setf(ios::fixed);
    cout.precision(2);

    benchmark_two_sum();
    benchmark_two_mult_FMA();

    benchmark_vec_sum(output_file);
    benchmark_vec_sum_err_branch(output_file);
    benchmark_vec_sum_err(output_file);

    benchmark_renormalization(output_file);
    benchmark_multiplication(output_file);

    output_file.close();
}
