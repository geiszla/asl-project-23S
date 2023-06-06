"""Module to plot benchmark results."""
import os
from dataclasses import dataclass
from os import path
from pathlib import Path
from typing import Literal, Union

import numpy
import pandas
from matplotlib import pyplot
from matplotlib.ticker import ScalarFormatter

from .algorithms import algorithms

BANDWIDTH = 18.7  # Bytes/cycle
PERFORMANCE_BOUND = 6  # flops/cycle
PERFORMANCE_BOUND_SIMD = 24  # flops/cycle

script_path = Path(__file__).parent.resolve()
results_path = path.join(script_path, "..", "results")


@dataclass
class Run:
    """Represents a run of the program."""

    label: str
    operational_intensity: float
    performance: float


def main():
    """Run script to generate performance plots."""
    os.makedirs(results_path, exist_ok=True)

    csv_data = pandas.read_csv(path.join(results_path, "benchmark.csv"), delimiter=";")

    runs: list[Run] = []
    reference_runs: list[Run] = []

    # Plot Input size vs Runtime columns for each Algorithm
    for algorithm in csv_data["Algorithm"].unique():
        optimization_data = csv_data[csv_data["Algorithm"] == algorithm]

        # Plot fast VecSum algorithms separately
        if algorithm == "VecSum":
            fast_data = optimization_data[
                optimization_data["Variant"].str.contains("Fast|Reference")
            ]
            optimization_data = optimization_data[
                # pylint: disable=singleton-comparison
                optimization_data["Variant"].str.contains("Fast")
                == False  # noqa: E712
            ]

            plot_performance("Performance", fast_data, "VecSumFast")
            plot_performance("Runtime", fast_data, "VecSumFast")

        plot_performance("Performance", optimization_data, algorithm)
        plot_performance("Runtime", optimization_data, algorithm)

        if algorithm not in algorithms:
            print(
                "Warning: could not find work and memory movement calculations"
                f" for {algorithm}"
            )
            continue

        # Calculate operational intensity
        # Assumes that the last variant of an algorithm is the reference, and the
        # last but one is the most optimized implementation
        variants = optimization_data["Variant"].unique()

        reference_data = optimization_data[optimization_data["Variant"] == variants[-1]]
        reference_runs.append(calculate_run_statistics(algorithm, reference_data))

        optimized_data = optimization_data[optimization_data["Variant"] == variants[-2]]
        runs.append(calculate_run_statistics(algorithm, optimized_data))

    # Create the roofline plot
    plot_roofline(runs)
    plot_roofline(runs, "Reference roofline model", "reference_roofline")


if __name__ == "__main__":
    main()


# Main functions


def plot_performance(
    performance_column: str, optimization_data: pandas.DataFrame, algorithm: str
):
    """Create performance plots using the given performance metrics."""
    pyplot.xscale("log", base=2)

    if performance_column == "Performance":
        pyplot.ylim([0, 2])

    pyplot.gca().get_xaxis().set_major_formatter(ScalarFormatter())

    variants = optimization_data["Variant"].unique()

    for variant in variants:
        variant_data = optimization_data[optimization_data["Variant"] == variant]
        input_sizes = variant_data["Input size"]
        performances = variant_data[performance_column]

        pyplot.plot(input_sizes, performances, label=variant)

    pyplot.legend()

    unit = "(flops/cycle)" if performance_column == "Performance" else "cycles"

    finalize_plot(
        x_label="Input size (n)",
        y_label=f"{performance_column} {unit}",
        title=algorithm,
        file_name=f"{performance_column.lower()}_vs_input_size_{algorithm}.png",
        is_show=True,
    )


def calculate_run_statistics(algorithm: str, optimization_data: pandas.DataFrame):
    """Calculate operational intensity and collect performance statistics."""
    sample_index = 5
    input_size = optimization_data["Input size"].iloc[sample_index]
    performance = optimization_data["Performance"].iloc[sample_index]

    flops = algorithms[algorithm].get_flops(input_size)
    memory_movements = algorithms[algorithm].get_memory_movements(input_size)
    operational_intensity = flops / memory_movements

    return Run(algorithm, operational_intensity, performance)


def plot_roofline(
    runs: list[Run], title: str = "Roofline model", filename: str = "roofline"
):
    """Draw the roofline plot and place the program runs on it."""
    # Draw performance bounds
    plot_horizontal_line(
        PERFORMANCE_BOUND, f"π bound without SIMD ({PERFORMANCE_BOUND})", "log"
    )
    plot_horizontal_line(
        PERFORMANCE_BOUND_SIMD, f"π bound with SIMD ({PERFORMANCE_BOUND_SIMD})", "log"
    )

    limit = PERFORMANCE_BOUND / BANDWIDTH
    plot_vertical_line(limit, f"π/β without SIMD ({limit:.2f})", "log")

    # Draw bandwidth bound
    x_coordinates = numpy.logspace(numpy.log10(0.001), numpy.log10(10), 100)
    y_coordinates = x_coordinates * BANDWIDTH

    pyplot.text(
        x_coordinates[10] * 2,
        y_coordinates[10],
        f"bound based on β ({BANDWIDTH} * I)",
        color="tab:blue",
    )

    pyplot.plot(x_coordinates, y_coordinates, color="tab:blue")
    pyplot.xlim([10e-4, 100])
    pyplot.xscale("log", base=2)
    pyplot.yscale("log", base=2)

    for run in runs:
        pyplot.plot(run.operational_intensity, run.performance, "o", label=run.label)

    pyplot.legend()

    finalize_plot(
        "Operational Intensity (flops/byte)",
        "Performance (flops/cycle)",
        title,
        f"{filename}.png",
        is_show=True,
    )


# Plotting helper functions


def plot_vertical_line(
    x_position: float, label: str, scale: Union[Literal["log"], Literal["linear"]]
):
    """Plot a vertical line at the given x position."""
    pyplot.axvline(x=x_position, color="tab:orange")  # L1 cache
    pyplot.text(
        x_position - (50 if scale == "linear" else x_position * 0.4),
        20,
        label,
        rotation=90,
        color="tab:orange",
    )


def plot_horizontal_line(
    y_position: float, label: str, scale: Union[Literal["log"], Literal["linear"]]
):
    """Plot a horizontal line at the given y position."""
    pyplot.axhline(y=y_position, color="tab:red")  # L1 cache
    pyplot.text(
        10e-4 * 1.5,
        y_position - (50 if scale == "linear" else y_position * 0.35),
        label,
        color="tab:red",
    )


def finalize_plot(
    x_label: float, y_label: float, title: str, file_name: str, is_show=False
):
    """Finalize the plot and save it to a file."""
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.title(title)
    pyplot.savefig(path.join(results_path, file_name))

    if is_show:
        pyplot.show()

    pyplot.clf()
