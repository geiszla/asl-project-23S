"""Module to plot benchmark results."""
import math
import os
from dataclasses import dataclass
from os import path
from pathlib import Path
from typing import Literal, Union

import numpy
import pandas
from matplotlib import pyplot
from matplotlib.ticker import ScalarFormatter

from .algorithms import memory_movement_calculations

BANDWIDTH = 18.7  # Bytes/cycle
PERFORMANCE_BOUND = 6  # flops/cycle
PERFORMANCE_BOUND_SIMD = 24  # flops/cycle

# Change it to, e.g., "Addition" or "Multiplication" to only show those graphs
ALGORITHM_FILTER = "???"

COLORS = ["gray", "purple", "blue", "orange", "red", "cyan", "black", "green"]
MARKERS = ["P", "D", "v", "^", "s", "x", "H", "o"]

script_path = Path(__file__).parent.resolve()
results_path = path.join(script_path, "..", "results")


@dataclass
class Run:
    """Represents a run of the program."""

    algorithm: str
    variant: str
    operational_intensity: float
    performance: float


def main():
    """Run script to generate performance plots."""
    os.makedirs(results_path, exist_ok=True)

    csv_data = pandas.read_csv(path.join(results_path, "benchmark.csv"), delimiter=";")

    runs: list[Run] = []

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

            runs.append(
                plot_optimizations(
                    fast_data,
                    algorithm,
                    filename="roofline_VecSumFast",
                )
            )

        plot_performance("Performance", optimization_data, algorithm)
        plot_performance("Runtime", optimization_data, algorithm)

        runs.append(plot_optimizations(optimization_data, algorithm))

        if algorithm not in memory_movement_calculations:
            print(
                "Warning: could not find work and memory movement calculations"
                f" for {algorithm}"
            )
            continue

    # Create the roofline plot
    plot_roofline(runs, is_connect=False)


# Main functions


def plot_performance(
    performance_column: str, optimization_data: pandas.DataFrame, algorithm: str
):
    """Create performance plots using the given performance metrics."""
    pyplot.xscale("log", base=2)

    if performance_column == "Performance":
        pyplot.ylim([0, 0.8 if max(optimization_data["Performance"]) < 0.8 else 2.1])

    pyplot.gca().get_xaxis().set_major_formatter(ScalarFormatter())

    variants = optimization_data["Variant"].unique()
    baseline = optimization_data[optimization_data["Variant"] == algorithm]

    for index, variant in enumerate(variants):
        variant_data = optimization_data[optimization_data["Variant"] == variant]

        is_reference = index == len(variants) - 1

        biggest_speedup = (
            max(baseline["Runtime"] / variant_data["Runtime"].values)
            if len(baseline) != 0
            else 1
        )

        label = variant if variant != algorithm else "Baseline"
        label = f"{algorithm} (CAMPARY)" if is_reference else label

        label += (
            f" ({biggest_speedup:.1f}x)"
            if not math.isnan(biggest_speedup) and biggest_speedup != 1
            else ""
        )

        pyplot.plot(
            variant_data["Input size"],
            variant_data[performance_column],
            label=label,
            color=COLORS[index] if not is_reference else "green",
            marker=MARKERS[index] if not is_reference else "o",
            linestyle="-" if not is_reference else "--",
        )

    unit = "[flops/cycle]" if performance_column == "Performance" else "cycles"

    handles, labels = pyplot.gca().get_legend_handles_labels()
    pyplot.legend(reversed(handles), reversed(labels))

    finalize_plot(
        x_label="Input size",
        y_label=f"{performance_column} {unit}",
        title="",
        file_name=f"{performance_column.lower()}_vs_input_size_{algorithm}.pdf",
        is_show=ALGORITHM_FILTER in algorithm if ALGORITHM_FILTER is not None else True,
    )


def plot_optimizations(
    optimization_data: pandas.DataFrame,
    algorithm: str,
    filename: str = None,
):
    """Calculate operational intensity and performance of the algorithm variants

    Assumes that the last variant of an algorithm is the reference, the first one is
    the least, and the last but one is the most optimized implementation
    """
    variants = optimization_data["Variant"].unique()

    variant_runs: list[Run] = []

    for variant in variants:
        variant_data = optimization_data[optimization_data["Variant"] == variant]

        sample_index = 6
        input_size = variant_data["Input size"].iloc[sample_index]
        performance = variant_data["Performance"].iloc[sample_index]
        flops = variant_data["Flop count"].iloc[sample_index]

        memory_movements = memory_movement_calculations[algorithm](input_size)
        operational_intensity = flops / memory_movements

        run = Run(algorithm, variant, operational_intensity, performance)
        variant_runs.append(run)

    plot_roofline(
        variant_runs,
        filename=f"roofline_{algorithm}" if filename is None else filename,
        is_show=ALGORITHM_FILTER in algorithm if ALGORITHM_FILTER is not None else True,
    )

    return variant_runs[-2]


def plot_roofline(
    runs: list[Run],
    is_connect=True,
    filename: str = "roofline",
    is_show=True,
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
    plot_vertical_line(limit)

    # Draw bandwidth bound
    x_coordinates = numpy.logspace(numpy.log10(0.001), numpy.log10(10), 100)
    y_coordinates = x_coordinates * BANDWIDTH

    pyplot.text(
        x_coordinates[10] * 1.8,
        y_coordinates[10],
        f"bound based on β ({BANDWIDTH} * I)",
        color="tab:blue",
    )

    pyplot.plot(x_coordinates, y_coordinates, color="tab:blue")
    pyplot.xlim([10e-4, 100])
    pyplot.xscale("log", base=2)
    pyplot.yscale("log", base=2)

    if is_connect:
        pyplot.plot(
            [run.operational_intensity for run in runs[:-1]],
            [run.performance for run in runs[:-1]],
            color="blue",
            marker="^",
            label=runs[0].algorithm,
        )

        pyplot.plot(
            runs[-1].operational_intensity,
            runs[-1].performance,
            color="green",
            marker="o",
            label=f"{runs[-1].algorithm} (CAMPARY)",
        )
    else:
        # Exclude non-fast variant of VecSum from the plots
        for index, run in enumerate(runs):
            if index != 1:
                pyplot.plot(
                    run.operational_intensity,
                    run.performance,
                    marker="^",
                    color=COLORS[index],
                    label=run.algorithm,
                )

    pyplot.legend()

    finalize_plot(
        "Operational Intensity [flops/byte]",
        "Performance [flops/cycle]",
        title="",
        file_name=f"{filename}.pdf",
        is_show=is_show,
    )


# Plotting helper functions


def plot_vertical_line(
    x_position: float,
    label: str = "",
    scale: Union[Literal["log"], Literal["linear"]] = "log",
):
    """Plot a vertical line at the given x position."""
    pyplot.axvline(x=x_position, color="tab:orange")  # L1 cache

    if label != "":
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
        2,
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
    pyplot.grid(axis="y")
    pyplot.tight_layout()
    pyplot.title(title)
    pyplot.savefig(path.join(results_path, file_name))

    if is_show:
        pyplot.show()

    pyplot.clf()


if __name__ == "__main__":
    main()
