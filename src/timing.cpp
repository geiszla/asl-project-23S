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

double rdtsc(comp_func f, quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N])
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
            f(x, y, A);
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
            f(x, y, A);
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

double query_performance_counter(comp_func f, quaternion_t x[N], quaternion_t y[N],
                                 quaternion_t A[N][N])
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
            f(x, y, A);
        }
        QueryPerformanceCounter(&end);

        cycles = (double)(end.QuadPart - start.QuadPart);
        multiplier = (CYCLES_REQUIRED) / (cycles * (FREQUENCY / frequency.QuadPart));
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
            f(x, y, A);
        }
        QueryPerformanceCounter(&end);

        cycles = (double)(end.QuadPart - start.QuadPart) / num_runs;
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= MEASUREMENT_COUNT;

    cycles = total_cycles;
    return cycles;
}

#endif

double measure_runtime(comp_func f, quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N])
{
#ifdef rdtsc
    return rdtsc(f, x, y, A);
#else
    return query_performance_counter(f, x, y, A);
#endif
}
