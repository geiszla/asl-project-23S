"""Module for calculations for different algorithms."""
# pylint: disable=unnecessary-lambda

from typing import Protocol

TWO_SUM_FLOPS = 6
FAST_TWO_SUM_FLOPS = 3
TWO_MULT_FMA_FLOPS = 3


# pylint: disable=too-few-public-methods
class FlopCalculation(Protocol):
    """Type of the flop and memory movement getter functions"""

    def __call__(self, *args: int) -> int:
        ...


# Algorithms

memory_movement_calculations: dict[str, FlopCalculation] = {
    "VecSum": lambda input_size: 8 * 3 * input_size,
    "VecSumErrBranch": lambda input_size: 8 * (input_size + 2 * input_size),
    "VecSumErr": lambda input_size: 8 * 3 * input_size,
    "Renormalization": lambda input_size: 8 * (input_size + 2 * input_size),
    "Addition": lambda input_size: 8 * (input_size + input_size + input_size),
    "Multiplication": lambda input_size: 8 * 3 * input_size,
    "FourMultiplications": lambda input_size: 4 * 8 * 3 * input_size,
    "Multiplication2": lambda input_size: 8 * 3 * input_size,
}
