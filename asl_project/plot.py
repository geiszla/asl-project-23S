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

    data_points: list[Run] = []

    # Plot Input size vs Runtime columns for each Algorithm
    for algorithm in csv_data["Algorithm"].unique():
        optimization_data = csv_data[csv_data["Algorithm"] == algorithm]
        input_sizes = optimization_data["Input size"]
        performances = optimization_data["Performance"]

        pyplot.xscale("log", base=2)
        pyplot.ylim(bottom=0)
        pyplot.gca().get_xaxis().set_major_formatter(ScalarFormatter())

        pyplot.plot(input_sizes, performances, label=algorithm)

        finalize_plot(
            x_label="Input size (n)",
            y_label="Performance (flops/cycle)",
            title=algorithm,
            file_name=f"performance_vs_input_size_{algorithm}.png",
        )

        if algorithm not in algorithms:
            print(
                "Warning: could not find work and memory movement calculations"
                f" for {algorithm}"
            )
            continue

        # Calculate operational intensity
        sample_index = 5
        input_size = input_sizes.iloc[sample_index]
        performance = performances.iloc[sample_index]

        flops = algorithms[algorithm].get_flops(input_size)
        memory_movements = algorithms[algorithm].get_memory_movements(input_size)
        operational_intensity = flops / memory_movements

        data_points.append(Run(algorithm, operational_intensity, performance))

    # Create the roofline plot
    plot_roofline(data_points)


def plot_roofline(runs: list[Run]):
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
        "Roofline model",
        "roofline.png",
        is_show=True,
    )


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


if __name__ == "__main__":
    main()
