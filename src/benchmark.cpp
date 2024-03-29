#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstring>

extern "C"
{
#include "reference.h"
#include "basefunctions.c"
#include "vec_sum_optimizations.c"
#include "vec_sum_err_optimizations.c"
#include "vec_sum_err_branch_optimizations.c"
#include "renormalization_optimizations.c"
#include "addition_optimizations.c"
#include "multiplication_optimizations.c"
#include "fourmult.c"
#include "mult2_optimizations.c"
}

#include "measure.cpp"
#include "reference_cpp.h"

// Compile with `g++ -std=c++17 ./benchmark.cpp` (plus additional optimization flags)

#define OUTPUT_PATH "../results"

// Flop counts

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

unsigned int get_optimized_renormalization_flops(int input_expansion_length, int output_expansion_length)
{
    int vec_sum_flops = get_vec_sum_flops(input_expansion_length);
    int vec_sum_err_branch_flops = (input_expansion_length - 2) * TWO_SUM_FLOPS;
    int vec_sum_err_flops = ((output_expansion_length - 2) * (output_expansion_length + 1) / 2.0) * TWO_SUM_FLOPS;

    return vec_sum_flops + vec_sum_err_branch_flops + vec_sum_err_flops;
}

unsigned int get_addition_flops(int a_length, int b_length, int result_length)
{
    return get_renormalization_flops(a_length + b_length, result_length);
}

unsigned int get_optimized_addition_flops(int a_length, int b_length, int result_length)
{
    return get_optimized_renormalization_flops(a_length + b_length, result_length);
}

unsigned int get_multiplication_flops(int expansion_length)
{
    int flops = TWO_MULT_FMA_FLOPS;

    for (int n = 0; n < expansion_length; n++)
    {
        for (int i = 0; i <= n; i++)
        {
            flops += TWO_MULT_FMA_FLOPS;
        }

        flops += get_vec_sum_flops(n * n + n + 1);
    }

    flops += (expansion_length - 1) * 2 + expansion_length * expansion_length;
    flops += get_renormalization_flops(expansion_length + 1, expansion_length);

    return flops;
}

unsigned int get_optimized_multiplication_flops(int expansion_length)
{
    int flops = TWO_MULT_FMA_FLOPS;

    for (int n = 0; n < expansion_length; n++)
    {
        for (int i = 0; i <= n; i++)
        {
            flops += TWO_MULT_FMA_FLOPS;
        }

        flops += get_vec_sum_flops(n * n + n + 1);
    }

    flops += (expansion_length - 1) * 2 + expansion_length * expansion_length;
    flops += get_optimized_renormalization_flops(expansion_length + 1, expansion_length);

    return flops;
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
    int vec_sum_err_flops = ((output_expansion_length - 2) * (output_expansion_length + 1) / 2.0) * FAST_TWO_SUM_FLOPS;

    return vec_sum_flops + vec_sum_err_branch_flops + vec_sum_err_flops;
}

unsigned int get_fast_addition_flops(int a_length, int b_length, int result_length)
{
    return get_fast_renormalization_flops(a_length + b_length, result_length);
}

unsigned int get_fast_multiplication_flops(int expansion_length)
{
    int flops = TWO_MULT_FMA_FLOPS;

    for (int n = 0; n < expansion_length; n++)
    {
        for (int i = 0; i <= n; i++)
        {
            flops += TWO_MULT_FMA_FLOPS;
        }

        flops += get_fast_vec_sum_flops(n * n + n + 1);
    }

    flops += (expansion_length - 1) * 2 + expansion_length * expansion_length;
    flops += get_fast_renormalization_flops(expansion_length + 1, expansion_length);

    return flops;
}

template <int Term_count>
unsigned int get_multiplication_reference_flops()
{
    int multiplication_flops = calculate_multiplication_flops<Term_count, Term_count, Term_count>();
    int renormalization_flops = get_fast_renormalization_flops(Term_count + 1, Term_count);

    return multiplication_flops + renormalization_flops;
}

unsigned int get_truncatedMul_flops(int expansion_length)
{
    int FPmul_flops = 1;
    int FPadd_flops = 1;

    int frexp_flops = 1;
    int ldexp_flops = 1;
    int floor_flops = 1;
    int division_flops = 1;

    int len_B = floor(expansion_length * 53 / 45 + 2);
    return division_flops + 1 + expansion_length * 2 * frexp_flops + 2 * ldexp_flops + (len_B + 1) * FPmul_flops + pow(expansion_length, 2) / 2 * (2 * FAST_TWO_SUM_FLOPS + 2 * FPadd_flops + division_flops) + (expansion_length + 1) * (1 * FPmul_flops + 2 * FAST_TWO_SUM_FLOPS + 1 * FPadd_flops + division_flops) + len_B * FPadd_flops + len_B * FAST_TWO_SUM_FLOPS;
}

// Helpers

void write_results(string algorithm_name, string variant_name, double runtime, double performance,
                   int flop_count, int input_size, ostream &output_file)
{
    output_file << algorithm_name << ";" << variant_name << ";" << input_size << ";" << runtime
                << ";" << performance << ";" << flop_count << endl;

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

void benchmark_vec_sum(
    ofstream &output_file,
    void (*implementation)(double *, double *, int) = vecSum,
    string variant_name = "VecSum",
    unsigned int (*get_flops)(int) = get_vec_sum_flops)
{
    cout << endl
         << variant_name << ":" << endl;

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
    {
        double runtime = measure_vec_sum(implementation, term_count);
        double flop_count = get_flops(term_count);
        double performance = flop_count / runtime;

        write_results("VecSum", variant_name, runtime, performance, flop_count, term_count,
                      output_file);
    }
}

void benchmark_vec_sum_err_branch(
    ofstream &output_file,
    void (*implementation)(double *, int, int, double *) = vecSumErrBranch,
    string variant_name = "VecSumErrBranch",
    unsigned int (*get_flops)(int) = get_vec_sum_err_branch_flops)
{
    cout << endl
         << variant_name << ":" << endl;

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
    {
        double runtime = measure_vec_sum_err_branch(implementation, term_count);
        double flop_count = get_flops(term_count);
        double performance = flop_count / runtime;

        write_results("VecSumErrBranch", variant_name, runtime, performance, flop_count, term_count,
                      output_file);
    }
}

void benchmark_vec_sum_err(
    ofstream &output_file,
    void (*implementation)(double *, int, double *) = vecSumErr,
    string variant_name = "VecSumErr",
    unsigned int (*get_flops)(int) = get_vec_sum_err_flops)
{
    cout << endl
         << variant_name << ":" << endl;

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
    {
        double runtime = measure_vec_sum_err(implementation, term_count);
        double flop_count = get_flops(term_count);
        double performance = flop_count / runtime;

        write_results("VecSumErr", variant_name, runtime, performance, flop_count, term_count,
                      output_file);
    }
}

void benchmark_renormalization(
    ofstream &output_file,
    void (*implementation)(double *, int, double *, int) = renormalizationalgorithm,
    string variant_name = "Renormalization",
    unsigned int (*get_flops)(int, int) = get_renormalization_flops)
{
    cout << endl
         << variant_name << ":" << endl;

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
    {
        double runtime = measure_renormalization(implementation, term_count);
        double flop_count = get_flops(term_count, term_count);
        double performance = flop_count / runtime;

        write_results("Renormalization", variant_name, runtime, performance, flop_count, term_count,
                      output_file);
    }
}

void benchmark_addition(
    ofstream &output_file,
    void (*implementation)(double *, double *, double *, int, int, int) = addition,
    string variant_name = "Addition",
    unsigned int (*get_flops)(int, int, int) = get_addition_flops)
{
    cout << endl
         << variant_name << ":" << endl;

    for (int term_count = 3; term_count <= 742; term_count *= 1.4)
    {
        double runtime = measure_addition(implementation, term_count);
        double flop_count = get_flops(term_count, term_count, term_count);
        double performance = flop_count / runtime;

        write_results("Addition", variant_name, runtime, performance, flop_count, term_count,
                      output_file);
    }
}

void benchmark_multiplication(
    ofstream &output_file,
    void (*implementation)(double *, double *, double *, int, int, int) = multiplication,
    string variant_name = "Multiplication",
    unsigned int (*get_flops)(int) = get_multiplication_flops)
{
    cout << endl
         << variant_name << ":" << endl;

    for (int term_count = 3; term_count <= 50; term_count += 2)
    {
        double runtime = measure_multiplication(implementation, term_count);
        double flop_count = get_flops(term_count);
        double performance = flop_count / runtime;

        write_results("Multiplication", variant_name, runtime, performance, flop_count, term_count,
                      output_file);
    }
}

void benchmark_fourmultiplications(
    ofstream &output_file,
    void (*implementation)(double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           const int, const int, const int),
    string variant_name = "FourMultiplications",
    unsigned int (*get_flops)(int) = get_fast_multiplication_flops)
{
    cout << endl
         << variant_name << ":" << endl;

    for (int term_count = 3; term_count <= 32; term_count += 2)
    {
        double runtime = measure_fourmultiplications(implementation, term_count);
        double flop_count = 4 * get_flops(term_count);
        double performance = flop_count / runtime;

        write_results("Multiplication", variant_name, runtime, performance, flop_count, term_count,
                      output_file);
    }
}

void benchmark_multiplication2(
    ofstream &output_file,
    void (*implementation)(double *, double *, double *, int, int, int) = mult2_fast,
    string variant_name = "Multiplication2",
    unsigned int (*get_flops)(int) = get_truncatedMul_flops)
{
    cout << endl
         << variant_name << ":" << endl;

    for (int term_count = 1; term_count < 30; term_count += 4)
    {
        double runtime = measure_multiplication2(implementation, term_count);
        double flop_count = get_flops(term_count);
        double performance = flop_count / runtime;

        write_results("Multiplication2", variant_name, runtime, performance, flop_count, term_count,
                      output_file);
    }
}

// References

template <int Term_count>
void measure_vec_sum_reference(ofstream &output_file)
{
    double runtime = measure_vec_sum(fast_VecSum<Term_count>, Term_count);
    double flop_count = get_fast_vec_sum_flops(Term_count);
    double performance = flop_count / runtime;

    write_results("VecSum", "VecSumReference", runtime, performance, flop_count, Term_count,
                  output_file);
}

void benchmark_vec_sum_reference(ofstream &output_file)
{
    cout << endl
         << "VecSumReference: " << endl;

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

template <int Term_count>
void measure_renormalization_reference(ofstream &output_file)
{
    double runtime = measure_renormalization(fast_renorm3L<Term_count, Term_count>, Term_count);
    double flop_count = get_fast_renormalization_flops(Term_count, Term_count);
    double performance = flop_count / runtime;

    write_results("Renormalization", "RenormalizationReference", runtime, performance, flop_count,
                  Term_count, output_file);
}

void benchmark_renormalization_reference(ofstream &output_file)
{
    cout << endl
         << "RenormalizationReference: " << endl;

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
    double flop_count = get_fast_addition_flops(Term_count, Term_count, Term_count);
    double performance = flop_count / runtime;

    write_results("Addition", "AdditionReference", runtime, performance, flop_count, Term_count,
                  output_file);
}

void benchmark_addition_reference(ofstream &output_file)
{
    cout << endl
         << "AdditionReference: " << endl;

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
    double flop_count = get_multiplication_reference_flops<Term_count>();
    double performance = flop_count / runtime;

    write_results("Multiplication", "MultiplicationReference", runtime, performance, flop_count,
                  Term_count, output_file);
}

void benchmark_multiplication_reference(ofstream &output_file)
{
    cout << endl
         << "MultiplicationReference: " << endl;

    measure_multiplication_reference<3>(output_file);
    measure_multiplication_reference<5>(output_file);
    measure_multiplication_reference<7>(output_file);
    measure_multiplication_reference<9>(output_file);
    measure_multiplication_reference<11>(output_file);
    measure_multiplication_reference<13>(output_file);
    measure_multiplication_reference<15>(output_file);
    measure_multiplication_reference<17>(output_file);
    measure_multiplication_reference<19>(output_file);
    measure_multiplication_reference<21>(output_file);
    measure_multiplication_reference<23>(output_file);
    measure_multiplication_reference<25>(output_file);
    measure_multiplication_reference<27>(output_file);
    measure_multiplication_reference<29>(output_file);
    measure_multiplication_reference<31>(output_file);
    measure_multiplication_reference<33>(output_file);
    measure_multiplication_reference<35>(output_file);
    measure_multiplication_reference<37>(output_file);
    measure_multiplication_reference<39>(output_file);
    measure_multiplication_reference<41>(output_file);
    measure_multiplication_reference<43>(output_file);
    measure_multiplication_reference<45>(output_file);
    measure_multiplication_reference<47>(output_file);
    measure_multiplication_reference<49>(output_file);
}

template <int K, int L, int M>
void measure_mult2_reference(ofstream &output_file)
{
    double runtime = measure_multiplication2(truncatedMul<K, L, M>, K);

    double flop_count = get_truncatedMul_flops(K);
    double performance = flop_count / runtime;

    write_results("Multiplication2", "Multiplication2Reference", runtime, performance, flop_count,
                  K, output_file);
}

void benchmark_multiplication2_reference(ofstream &output_file)
{
    cout << endl
         << "Multiplication2 reference: " << endl;

    measure_mult2_reference<1, 1, 1>(output_file);
    measure_mult2_reference<5, 5, 5>(output_file);
    measure_mult2_reference<9, 9, 9>(output_file);
    measure_mult2_reference<13, 13, 13>(output_file);
    measure_mult2_reference<17, 17, 17>(output_file);
    measure_mult2_reference<21, 21, 21>(output_file);
    measure_mult2_reference<25, 25, 25>(output_file);
    measure_mult2_reference<29, 29, 29>(output_file);
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

    output_file << "Algorithm;Variant;Input size;Runtime;Performance;Flop count" << endl;

    cout.setf(ios::fixed);
    cout.precision(3);

    // TwoSum
    benchmark_two_sum();

    // TwoMultFMA
    benchmark_two_mult_FMA();

    // VecSum
    benchmark_vec_sum(output_file);
    benchmark_vec_sum(output_file, vecSum2, "VecSum2");
    benchmark_vec_sum(output_file, vecSum3, "VecSum3");
    benchmark_vec_sum(output_file, vecSum4, "VecSum4");
    benchmark_vec_sum(output_file, vecSum5, "VecSum5");
    benchmark_vec_sum(output_file, vecSum6, "VecSum6");
    benchmark_vec_sum(output_file, vecSum3_fast, "VecSum3Fast", get_fast_vec_sum_flops);
    benchmark_vec_sum(output_file, vecSum4_fast, "VecSum4Fast", get_fast_vec_sum_flops);
    benchmark_vec_sum(output_file, vecSum5_fast, "VecSum5Fast", get_fast_vec_sum_flops);
    benchmark_vec_sum(output_file, vecSum6_fast, "VecSum6Fast", get_fast_vec_sum_flops);
    benchmark_vec_sum_reference(output_file);

    // VecSumErrBranch
    benchmark_vec_sum_err_branch(output_file);
    benchmark_vec_sum_err_branch(output_file, vecSumErrBranch2, "VecSumErrBranch2");
    benchmark_vec_sum_err_branch(output_file, vecSumErrBranch_fast, "VecSumErrBranchFast",
                                 get_fast_vec_sum_err_branch_flops);
    benchmark_vec_sum_err_branch(output_file, fast_VecSumErrBranch, "VecSumErrBranchReference",
                                 get_fast_vec_sum_err_branch_flops);

    // VecSumErr
    benchmark_vec_sum_err(output_file);
    benchmark_vec_sum_err(output_file, vecSumErr2, "VecSumErr2");
    benchmark_vec_sum_err(output_file, vecSumErr_fast, "VecSumErrFast",
                          get_fast_vec_sum_err_flops);
    benchmark_vec_sum_err(output_file, fast_VecSumErr, "VecSumErrReference",
                          get_fast_vec_sum_err_flops);

    // Renormalization
    benchmark_renormalization(output_file);
    benchmark_renormalization(output_file, renormalization2, "Renormalization2");
    benchmark_renormalization(output_file, renormalization3, "Renormalization3");
    benchmark_renormalization(output_file, renormalization4, "Renormalization4",
                              get_optimized_renormalization_flops);
    benchmark_renormalization(output_file, renormalization_fast, "RenormalizationFast",
                              get_fast_renormalization_flops);
    benchmark_renormalization(output_file, renormalization_fast2, "RenormalizationFast2",
                              get_fast_renormalization_flops);
    benchmark_renormalization_reference(output_file);

    // Addition
    benchmark_addition(output_file);
    benchmark_addition(output_file, addition2, "Addition2");
    benchmark_addition(output_file, addition3, "Addition3");
    benchmark_addition(output_file, addition4, "Addition4", get_optimized_addition_flops);
    benchmark_addition(output_file, addition_fast, "AdditionFast", get_fast_addition_flops);
    benchmark_addition(output_file, addition_fast2, "AdditionFast2", get_fast_addition_flops);
    benchmark_addition_reference(output_file);

    // Multiplication
    benchmark_multiplication(output_file);
    benchmark_multiplication(output_file, multiplication2, "Multiplication2",
                             get_optimized_multiplication_flops);
    benchmark_multiplication(output_file, multiplication3, "Multiplication3",
                             get_optimized_multiplication_flops);

    // benchmark_fourmultiplications(output_file, fourtimesmultiplicationversion0,
    //                               "FourMultiplications0");
    // benchmark_fourmultiplications(output_file, fourtimesmultiplicationversion1,
    //                               "FourMultiplications1");
    // benchmark_fourmultiplications(output_file, fourtimesmultiplicationversion3,
    //                               "FourMultiplications3");

    benchmark_multiplication(output_file, multiplication_fast, "MultiplicationFast",
                             get_fast_multiplication_flops);
    benchmark_multiplication_reference(output_file);

    // Multiplication2
    benchmark_multiplication2(output_file);
    benchmark_multiplication2(output_file, mult2_0_fast, "Multiplication2_0");
    benchmark_multiplication2(output_file, mult2_1_fast, "Multiplication2_1");
    benchmark_multiplication2(output_file, mult2_2_fast, "Multiplication2_2");
    benchmark_multiplication2(output_file, mult2_3_fast, "Multiplication2_3");
    benchmark_multiplication2_reference(output_file);

    output_file.close();
}
