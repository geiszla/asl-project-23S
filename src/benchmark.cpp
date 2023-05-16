#include <filesystem>
#include <fstream>
#include <iostream>

#include "measure.cpp"

// Compile with `g++ -std=c++17 ./benchmark.cpp` (plus additional optimization flags)
// TODO: calculate median of measurements, generalize benchmark logic

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

unsigned int get_multiplication2_flops(int expansion_length)
{ // return Flp count from paper 2, Proposition 3.7
    return 13 / 2 * pow(expansion_length, 2) + 33 / 2 * expansion_length + 6 * (floor(expansion_length * 53 / 45) + 2) + 55;
}

// Helpers

void write_results(string algorithm_name, double runtime, double performance, int input_size,
                   ostream &output_file)
{
    output_file << algorithm_name << ";" << input_size << ";" << runtime
                << ";" << performance << endl;

    cout << "Input size: " << input_size
         << "     runtime: " << runtime << " cycles, performance: "
         << performance << " flops/cycles" << endl;
}

// Benchmarking functions

void benchmark_two_sum()
{
    double runtime = measure_two_sum();

    cout << "2Sum\t\truntime: " << runtime << " cycles, performance: "
         << (TWO_SUM_FLOPS / runtime) << " flops/cycles" << endl;
}

void benchmark_two_mult_FMA()
{
    double runtime = measure_two_mult_FMA();

    cout << "2MultFMA\truntime: " << runtime << " cycles, performance: "
         << (TWO_MULT_FMA_FLOPS / runtime) << " flops/cycles" << endl;
}

void benchmark_vec_sum(ofstream &output_file)
{
    cout << endl
         << "VecSum: " << endl;

    for (int term_count = 3; term_count < 999; (int)(term_count *= 1.5))
    {
        double runtime = measure_vec_sum(vecSum, term_count);
        double performance = get_vec_sum_flops(term_count) / runtime;

        write_results("VecSum", runtime, performance, term_count, output_file);
    }
}

void benchmark_vec_sum_err_branch(ofstream &output_file)
{
    cout << endl
         << "VecSumErrBranch: " << endl;

    for (int term_count = 3; term_count < 999; (int)(term_count *= 1.4))
    {
        double runtime = measure_vec_sum_err_branch(vecSumErrBranch, term_count);
        double performance = get_vec_sum_err_branch_flops(term_count) / runtime;

        write_results("VecSumErrBranch", runtime, performance, term_count, output_file);
    }
}

void benchmark_vec_sum_err(ofstream &output_file)
{
    cout << endl
         << "VecSumErr: " << endl;

    for (int term_count = 3; term_count < 999; (int)(term_count *= 1.4))
    {
        double runtime = measure_vec_sum_err(vecSumErr, term_count);
        double performance = get_vec_sum_err_flops(term_count) / runtime;

        write_results("VecSumErr", runtime, performance, term_count, output_file);
    }
}

void benchmark_renormalization(ofstream &output_file)
{
    cout << endl
         << "Renormalization: " << endl;

    for (int term_count = 3; term_count < 999; (int)(term_count *= 1.4))
    {
        double runtime = measure_renormalization(renormalizationalgorithm, term_count);
        double performance = get_renormalization_flops(term_count, term_count) / runtime;

        write_results("Renormalization", runtime, performance, term_count, output_file);
    }
}

void benchmark_addition(ofstream &output_file)
{
    cout << endl
         << "Addition: " << endl;

    for (int term_count = 3; term_count < 999; (int)(term_count *= 1.4))
    {
        double runtime = measure_addition(addition, term_count);
        double performance = get_addition_flops(term_count, term_count, term_count) / runtime;

        write_results("Addition", runtime, performance, term_count, output_file);
    }
}

void benchmark_multiplication(ofstream &output_file)
{
    cout << endl
         << "Multiplication: " << endl;

    for (int term_count = 3; term_count < 999; (int)(term_count *= 1.4))
    {
        double runtime = measure_multiplication(multiplication, term_count);
        double performance = get_multiplication_flops(term_count) / runtime;

        write_results("Multiplication", runtime, performance, term_count, output_file);
    }
}

void benchmark_multiplication2(ofstream &output_file)
{
    cout << endl
         << "Multiplication2: " << endl;

    for (int term_count = 1; term_count < 129; term_count *= 2)
    {
        double runtime = measure_multiplication2(mult2, term_count);
        double performance = get_multiplication2_flops(term_count) / runtime;

        write_results("Multiplication2", runtime, performance, term_count, output_file);
    }
}

void benchmark_multiplication2_0(ofstream &output_file)
{
    cout << endl
         << "Multiplication2_0: " << endl;

    for (int term_count = 1; term_count < 129; term_count *= 2)
    {
        double runtime = measure_multiplication2(mult2_0, term_count);
        double performance = get_multiplication2_flops(term_count) / runtime;

        write_results("Multiplication2_0", runtime, performance, term_count, output_file);
    }
}

void benchmark_multiplication2_1(ofstream &output_file)
{
    cout << endl
         << "Multiplication2_1: " << endl;

    for (int term_count = 1; term_count < 129; term_count *= 2)
    {
        double runtime = measure_multiplication2(mult2_1, term_count);
        double performance = get_multiplication2_flops(term_count) / runtime;

        write_results("Multiplication2_1", runtime, performance, term_count, output_file);
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
    cout.precision(3);

    benchmark_two_sum();
    benchmark_two_mult_FMA();

    benchmark_vec_sum(output_file);
    benchmark_vec_sum_err_branch(output_file);
    benchmark_vec_sum_err(output_file);

    // benchmark_renormalization(output_file);
    // benchmark_addition(output_file);
    // benchmark_multiplication(output_file);

    // benchmark_multiplication2(output_file);
    // benchmark_multiplication2_0(output_file);
    // benchmark_multiplication2_1(output_file);

    output_file.close();
}
