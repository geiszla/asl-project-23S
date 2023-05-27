#ifdef _WIN32
#include <windows.h> // Include if under windows
#endif

#include <algorithm>
#include <chrono>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef WIN32
#include <sys/time.h>
#endif

#ifdef __x86_64__
#include "tsc_x86.h"
#endif

using namespace std;

#define CALIBRATE
#define MEASUREMENT_COUNT 50
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2.9e9 // Change this to your base clock speed

/*
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__

template <typename T_Function, class... T_compute_arguments>
vector<double> rdtsc(T_Function f, T_compute_arguments... arguments)
{
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;

    myInt64 start, end;

    /*
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    do
    {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++)
        {
            f(arguments...);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    vector<double> cycles_list;

    // Actual performance measurements repeated MEASUREMENT_COUNT times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < MEASUREMENT_COUNT; j++)
    {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i)
        {
            f(arguments...);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;

        cycles_list.push_back(cycles);
    }

    return cycles_list;
}
#elif defined __MACH__

#define ctime_t std::chrono::high_resolution_clock::time_point

template <typename T_Function, class... T_compute_arguments>
vector<double> high_resolution_clock(T_Function f, T_compute_arguments... arguments)
{
    uint64_t frequency = 3.5e9;

    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;

    ctime_t start, end;

#ifdef CALIBRATE
    do
    {
        num_runs = num_runs * multiplier;

        start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < num_runs; ++i)
        {
            f(arguments...);
        }
        end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;
        cycles = elapsed.count() * frequency;
        multiplier = CYCLES_REQUIRED / cycles;
    } while (multiplier > 2);
#endif

    vector<double> cycles_list;

    // Actual performance measurements repeated MEASUREMENT_COUNT times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < MEASUREMENT_COUNT; j++)
    {
        start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < num_runs; ++i)
        {
            f(arguments...);
        }
        end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;
        cycles = elapsed.count() * frequency;

        cycles_list.push_back(cycles);
    }

    return cycles_list;
}

#else

template <typename T_Function, class... T_compute_arguments>
vector<double> query_performance_counter(T_Function f, T_compute_arguments... arguments)
{
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);

    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;

    LARGE_INTEGER start, end;

#ifdef CALIBRATE
    do
    {
        num_runs = num_runs * multiplier;

        QueryPerformanceCounter(&start);
        for (size_t i = 0; i < num_runs; ++i)
        {
            f(arguments...);
        }
        QueryPerformanceCounter(&end);

        cycles = (double)(end.QuadPart - start.QuadPart) * (FREQUENCY / frequency.QuadPart);
        multiplier = CYCLES_REQUIRED / cycles;
    } while (multiplier > 2);
#endif

    vector<double> cycles_list;

    // Actual performance measurements repeated MEASUREMENT_COUNT times.
    // We simply store all results and compute medians during post-processing.
    for (size_t j = 0; j < MEASUREMENT_COUNT; j++)
    {
        QueryPerformanceCounter(&start);
        for (size_t i = 0; i < num_runs; ++i)
        {
            f(arguments...);
        }
        QueryPerformanceCounter(&end);

        cycles = (double)(end.QuadPart - start.QuadPart) / num_runs;
        cycles *= (FREQUENCY / frequency.QuadPart);

        cycles_list.push_back(cycles);
    }

    return cycles_list;
}

#endif

template <typename T_compute_return, class... T_compute_arguments>
double measure_runtime(T_compute_return (*f)(T_compute_arguments...),
                       T_compute_arguments... arguments)
{
    vector<double> cycles_list;

#ifdef __x86_64__
    cycles_list = rdtsc(f, arguments...);
#elif defined __MACH__
    cycles_list = high_resolution_clock(f, arguments...);
#else
    cycles_list = query_performance_counter(f, arguments...);
#endif

    // Calculate median
    vector<double>::iterator middle_iterator = cycles_list.begin() + cycles_list.size() / 2;
    nth_element(cycles_list.begin(), middle_iterator, cycles_list.end());

    if (cycles_list.size() % 2 == 0)
    {
        return (*max_element(cycles_list.begin(), middle_iterator) + *middle_iterator) / 2;
    }

    return *middle_iterator;
}
