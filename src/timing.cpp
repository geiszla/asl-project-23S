#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main.h"

#ifndef WIN32
#include <sys/time.h>
#endif

#include "main.h"

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

template <typename T_compute_function, class... T_compute_arguments>
double query_performance_counter(T_compute_function f, T_compute_arguments... arguments)
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

    list<double> cyclesList;

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
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= MEASUREMENT_COUNT;

    cycles = total_cycles;
    return cycles;
}

#else

template <typename T_compute_function, class... T_compute_arguments>
double query_performance_counter(T_compute_function f, T_compute_arguments... arguments)
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

    list<double> cyclesList;

    // Actual performance measurements repeated MEASUREMENT_COUNT times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
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

        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= MEASUREMENT_COUNT;

    cycles = total_cycles;
    return cycles;
}

#endif

template <typename T_compute_return, class... T_compute_arguments>
double measure_runtime(compute_function<T_compute_return, T_compute_arguments...> f,
                                 T_compute_arguments... arguments)
{
#ifdef rdtsc
    return rdtsc(f, arguments...);
#else
    return query_performance_counter(f, arguments...);
#endif
}
