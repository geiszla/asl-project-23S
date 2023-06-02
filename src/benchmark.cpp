#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstring>

extern "C"
{
#include "reference.h"
#include "basefunctions.c"
}

#include "measure.cpp"
#include "reference_cpp.h"

// Compile with `g++ -std=c++17 ./benchmark.cpp` (plus additional optimization flags)
// TODO: fix multiplication flop count

#define OUTPUT_PATH "../results"

// Flop counts

#define TWO_SUM_FLOPS 6
#define FAST_TWO_SUM_FLOPS 3
#define TWO_MULT_FMA_FLOPS 2

// #define BENCHMARK_REFERENCES

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

// Reference flop counts

unsigned int get_fast_vec_sum_flops(int vector_length)
{
    return (vector_length - 1) * FAST_TWO_SUM_FLOPS;
}

unsigned int get_fast_vec_sum_err_branch_flops(int input_expansion_length)
{
    return (input_expansion_length - 1) * FAST_TWO_SUM_FLOPS;
}

unsigned int get_fast_vec_sum_err_flops(int expansion_length)
{
    return (expansion_length - 1) * FAST_TWO_SUM_FLOPS;
}

unsigned int get_fast_renormalization_flops(int input_expansion_length, int output_expansion_length)
{
    int vec_sum_flops = get_fast_vec_sum_flops(input_expansion_length);
    int vec_sum_err_branch_flops = (input_expansion_length - 2) * FAST_TWO_SUM_FLOPS;
    int vec_sum_err_flops = get_fast_vec_sum_err_flops(
        (output_expansion_length * (output_expansion_length + 1) / 2.0));

    return vec_sum_flops + vec_sum_err_branch_flops + vec_sum_err_flops;
}

unsigned int get_addition_reference_flops(int a_length, int b_length, int result_length)
{
    int input_expansion_length = a_length + b_length;

    int vec_sum_flops = get_vec_sum_flops(input_expansion_length);
    int vec_sum_err_branch_flops = (input_expansion_length - 2) * FAST_TWO_SUM_FLOPS;

    return vec_sum_flops + vec_sum_err_branch_flops;
}

unsigned int get_multiplication_reference_flops(int expansion_length)
{
    // TODO: fix flop count
    return expansion_length * (expansion_length + 1) / 2 * TWO_MULT_FMA_FLOPS +
           (expansion_length / 2) * get_vec_sum_flops(expansion_length) +
           (expansion_length - 1) * 2 +
           pow(expansion_length, 2) +
           get_renormalization_flops(expansion_length + 1, expansion_length);
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

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
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

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
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

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
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

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
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

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
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

    for (int term_count = 3; term_count <= 742; term_count += 2)
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

    for (int term_count = 1; term_count <= 29; term_count *= 2)
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

    for (int term_count = 1; term_count <= 29; term_count *= 2)
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

    for (int term_count = 1; term_count <= 29; term_count *= 2)
    {
        double runtime = measure_multiplication2(mult2_1, term_count);
        double performance = get_multiplication2_flops(term_count) / runtime;

        write_results("Multiplication2_1", runtime, performance, term_count, output_file);
    }
}

void benchmark_multiplication2_2(ofstream &output_file)
{
    cout << endl
         << "Multiplication2_2: " << endl;

    for (int term_count = 1; term_count <= 29; term_count *= 2)
    {
        double runtime = measure_multiplication2(mult2_1, term_count);
        double performance = get_multiplication2_flops(term_count) / runtime;

        write_results("Multiplication2_2", runtime, performance, term_count, output_file);
    }
}

void benchmark_multiplication2_3(ofstream &output_file)
{
    cout << endl
         << "Multiplication2_3: " << endl;

    for (int term_count = 1; term_count <= 29; term_count *= 2)
    {
        double runtime = measure_multiplication2(mult2_1, term_count);
        double performance = get_multiplication2_flops(term_count) / runtime;

        write_results("Multiplication2_3", runtime, performance, term_count, output_file);
    }
}

// References

template <int Term_count>
void measure_vec_sum_reference(ofstream &output_file)
{
    double runtime = measure_vec_sum(fast_VecSum<Term_count>, Term_count);
    double performance = get_fast_vec_sum_flops(Term_count) / runtime;

    write_results("VecSumReference", runtime, performance, Term_count, output_file);
}

void benchmark_vec_sum_reference(ofstream &output_file)
{
    cout << endl
         << "VecSum reference: " << endl;

    measure_vec_sum_reference<3>(output_file);
    measure_vec_sum_reference<4>(output_file);
    measure_vec_sum_reference<5>(output_file);
    measure_vec_sum_reference<7>(output_file);
    measure_vec_sum_reference<9>(output_file);
    measure_vec_sum_reference<12>(output_file);
    measure_vec_sum_reference<16>(output_file);
    measure_vec_sum_reference<22>(output_file);
    measure_vec_sum_reference<30>(output_file);
    measure_vec_sum_reference<42>(output_file);
    measure_vec_sum_reference<58>(output_file);
    measure_vec_sum_reference<81>(output_file);
    measure_vec_sum_reference<113>(output_file);
    measure_vec_sum_reference<158>(output_file);
    measure_vec_sum_reference<221>(output_file);
    measure_vec_sum_reference<309>(output_file);
    measure_vec_sum_reference<432>(output_file);
    measure_vec_sum_reference<604>(output_file);
}

void benchmark_vec_sum_err_branch_reference(ofstream &output_file)
{
    cout << endl
         << "VecSumErrBranch Reference: " << endl;

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
    {
        double runtime = measure_vec_sum_err_branch(fast_VecSumErrBranch, term_count);
        double performance = get_fast_vec_sum_err_branch_flops(term_count) / runtime;

        write_results("VecSumErrBranchReference", runtime, performance, term_count, output_file);
    }
}

void benchmark_vec_sum_err_reference(ofstream &output_file)
{
    cout << endl
         << "VecSumErr Reference: " << endl;

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
    {
        double runtime = measure_vec_sum_err(fast_VecSumErr, term_count);
        double performance = get_fast_vec_sum_err_flops(term_count) / runtime;

        write_results("VecSumErrReference", runtime, performance, term_count, output_file);
    }
}

template <int Term_count>
void measure_renormalization_reference(ofstream &output_file)
{
    double runtime = measure_renormalization(fast_renorm3L<Term_count, Term_count>, Term_count);
    double performance = get_fast_renormalization_flops(Term_count, Term_count) / runtime;

    write_results("RenormalizationReference", runtime, performance, Term_count, output_file);
}

void benchmark_renormalization_reference(ofstream &output_file)
{
    cout << endl
         << "Renormalization Reference: " << endl;

    measure_renormalization_reference<3>(output_file);
    measure_renormalization_reference<4>(output_file);
    measure_renormalization_reference<5>(output_file);
    measure_renormalization_reference<7>(output_file);
    measure_renormalization_reference<9>(output_file);
    measure_renormalization_reference<12>(output_file);
    measure_renormalization_reference<16>(output_file);
    measure_renormalization_reference<22>(output_file);
    measure_renormalization_reference<30>(output_file);
    measure_renormalization_reference<42>(output_file);
    measure_renormalization_reference<58>(output_file);
    measure_renormalization_reference<81>(output_file);
    measure_renormalization_reference<113>(output_file);
    measure_renormalization_reference<158>(output_file);
    measure_renormalization_reference<221>(output_file);
    measure_renormalization_reference<309>(output_file);
    measure_renormalization_reference<432>(output_file);
    measure_renormalization_reference<604>(output_file);
}

template <int Term_count>
void measure_addition_reference(ofstream &output_file)
{
    double runtime = measure_addition(certifiedAdd<Term_count, Term_count, Term_count>, Term_count);
    double performance = get_addition_reference_flops(Term_count, Term_count, Term_count) / runtime;

    write_results("AdditionReference", runtime, performance, Term_count, output_file);
}

void benchmark_addition_reference(ofstream &output_file)
{
    cout << endl
         << "Addition Reference: " << endl;

    measure_addition_reference<3>(output_file);
    measure_addition_reference<4>(output_file);
    measure_addition_reference<5>(output_file);
    measure_addition_reference<7>(output_file);
    measure_addition_reference<9>(output_file);
    measure_addition_reference<12>(output_file);
    measure_addition_reference<16>(output_file);
    measure_addition_reference<22>(output_file);
    measure_addition_reference<30>(output_file);
    measure_addition_reference<42>(output_file);
    measure_addition_reference<58>(output_file);
    measure_addition_reference<81>(output_file);
    measure_addition_reference<113>(output_file);
    measure_addition_reference<158>(output_file);
    measure_addition_reference<221>(output_file);
    measure_addition_reference<309>(output_file);
    measure_addition_reference<432>(output_file);
    measure_addition_reference<604>(output_file);
}

template <int Term_count>
void measure_multiplication_reference(ofstream &output_file)
{
    double runtime = measure_multiplication(baileyMul_renorm<Term_count, Term_count, Term_count>,
                                            Term_count);
    double performance = get_multiplication_reference_flops(Term_count) / runtime;

    write_results("MultiplicationReference", runtime, performance, Term_count, output_file);
}

void benchmark_multiplication_reference(ofstream &output_file)
{
    cout << endl
         << "Multiplication Reference: " << endl;

    measure_multiplication_reference<3>(output_file);
    measure_multiplication_reference<4>(output_file);
    measure_multiplication_reference<5>(output_file);
    measure_multiplication_reference<7>(output_file);
    measure_multiplication_reference<9>(output_file);
    measure_multiplication_reference<12>(output_file);
    measure_multiplication_reference<16>(output_file);
    measure_multiplication_reference<22>(output_file);
    measure_multiplication_reference<30>(output_file);
    measure_multiplication_reference<42>(output_file);
    measure_multiplication_reference<58>(output_file);
    measure_multiplication_reference<81>(output_file);
    measure_multiplication_reference<113>(output_file);
    measure_multiplication_reference<158>(output_file);
    measure_multiplication_reference<221>(output_file);
    measure_multiplication_reference<309>(output_file);
    measure_multiplication_reference<432>(output_file);
    measure_multiplication_reference<604>(output_file);
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

#if !defined(BENCHMARK_REFERENCES)
    benchmark_two_sum();
    benchmark_two_mult_FMA();

    benchmark_vec_sum(output_file);
    benchmark_vec_sum_err_branch(output_file);
    benchmark_vec_sum_err(output_file);

    benchmark_renormalization(output_file);
    benchmark_addition(output_file);
    benchmark_multiplication(output_file);

    benchmark_multiplication2(output_file);
    benchmark_multiplication2_0(output_file);
    benchmark_multiplication2_1(output_file);
    benchmark_multiplication2_2(output_file);
    benchmark_multiplication2_3(output_file);
#else
    benchmark_vec_sum_reference(output_file);
    benchmark_vec_sum_err_branch_reference(output_file);
    benchmark_vec_sum_err_reference(output_file);
    benchmark_renormalization_reference(output_file);
    benchmark_addition_reference(output_file);
    benchmark_multiplication_reference(output_file);
#endif

    output_file.close();
}
