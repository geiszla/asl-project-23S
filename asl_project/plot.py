"""Module to plot benchmark results."""
import os

import pandas
from matplotlib import pyplot

RESULT_DIRECTORY_NAME = "results"

CLOCK_SPEED = 2.9e9  # Hz


def main():
    """Run script to generate plots for homework 1."""
    os.makedirs(RESULT_DIRECTORY_NAME, exist_ok=True)

    csv_data = pandas.read_csv("results/benchmark.csv", delimiter=";")

    for algorithm in csv_data["Algorithm"].unique():
        # Plot Input vs Runtime columns for each Optimizations value
        optimization_data = csv_data[csv_data["Algorithm"] == algorithm]
        input_size = optimization_data["Input size"]
        performance = optimization_data["Performance"]

        pyplot.plot(
            input_size,
            performance,
            label=algorithm,
        )

        pyplot.legend()
        finalize_plot(
            x_label="Input size (n)",
            y_label="Performance (flops/cycle)",
            title="Performance vs Input size",
            file_name=f"performance_vs_input_size_{algorithm}.png",
        )


def finalize_plot(x_label, y_label, title, file_name):
    """Finalize the plot and save it to a file."""
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.title(title)
    pyplot.savefig(f"{RESULT_DIRECTORY_NAME}/{file_name}")
    pyplot.show()
    pyplot.clf()


if __name__ == "__main__":
    main()
