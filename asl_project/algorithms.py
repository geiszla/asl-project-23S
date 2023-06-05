"""Module for calculations for different algorithms."""
# pylint: disable=unnecessary-lambda

from dataclasses import dataclass
from typing import Protocol

from math import floor

TWO_SUM_FLOPS = 6
FAST_TWO_SUM_FLOPS = 3
TWO_MULT_FMA_FLOPS = 2


# pylint: disable=too-few-public-methods
class Getter(Protocol):
    """Type of the flop and memory movement getter functions"""

    def __call__(self, *args: int) -> int:
        ...


@dataclass
class Algorithm:
    """Represents and algorithm."""

    get_flops: Getter
    get_memory_movements: Getter


# Algorithms

algorithms = {
    "VecSum": Algorithm(
        get_flops=lambda input_size: (input_size - 1) * TWO_SUM_FLOPS,
        get_memory_movements=lambda input_size: 8 * 3 * input_size,
    ),
    "VecSumErrBranch": Algorithm(
        get_flops=lambda input_size: (input_size - 1) * TWO_SUM_FLOPS,
        get_memory_movements=lambda input_size: 8 * (input_size + 2 * input_size),
    ),
    "VecSumErr": Algorithm(
        get_flops=lambda input_size: (input_size - 1) * TWO_SUM_FLOPS,
        get_memory_movements=lambda input_size: 8 * 3 * input_size,
    ),
    "Renormalization": Algorithm(
        get_flops=lambda input_size: get_renormalization_flops(input_size, input_size),
        get_memory_movements=lambda input_size: 8 * (input_size + 2 * input_size),
    ),
    "Addition": Algorithm(
        get_flops=lambda input_size: get_renormalization_flops(
            input_size + input_size, input_size
        ),
        get_memory_movements=lambda input_size: 8
        * (input_size + input_size + input_size),
    ),
    "Multiplication": Algorithm(
        get_flops=lambda input_size: get_multiplication_flops(input_size),
        get_memory_movements=lambda input_size: 8 * 3 * input_size,
    ),
        "Multiplication2": Algorithm(
        get_flops=lambda input_size: get_multiplication2_flops(input_size),
        get_memory_movements=lambda input_size: 8 * 3 * input_size,
    ),
    # "VecSumReference": Algorithm(
    #     get_flops=lambda vector_length: (vector_length - 1) * FAST_TWO_SUM_FLOPS,
    #     get_memory_movements=lambda vector_length: 8 * 2 * vector_length,
    # ),
    # "VecSumErrBranchReference": Algorithm(
    #     get_flops=lambda input_length: (input_length - 1) * FAST_TWO_SUM_FLOPS,
    #     get_memory_movements=lambda input_length, output_length: 8
    #     * (input_length + output_length),
    # ),
    # "VecSumErrReference": Algorithm(
    #     get_flops=lambda expansion_length: (expansion_length - 1) * FAST_TWO_SUM_FLOPS,
    #     get_memory_movements=lambda expansion_length: 8 * 2 * expansion_length,
    # ),
    # "RenormalizationReference": Algorithm(
    #     get_flops=lambda input_length, output_length: get_fast_renormalization_flops(
    #         input_length, output_length
    #     ),
    #     get_memory_movements=lambda expansion_length: 8 * 2 * expansion_length,
    # ),
    # "AdditionReference": Algorithm(
    #     get_flops=lambda a_length, b_length: get_addition_reference_flops(
    #         a_length, b_length
    #     ),
    #     get_memory_movements=lambda expansion_length: 8 * 2 * expansion_length,
    # ),
    # "MultiplicationReference": Algorithm(
    #     get_flops=lambda expansion_length: get_multiplication_reference_flops(
    #         expansion_length
    #     ),
    #     get_memory_movements=lambda expansion_length: 8 * 2 * expansion_length,
    # ),
}


# Helpers


def get_renormalization_flops(input_length: int, output_length: int):
    """Calculate the number of floating point operations done by Renormalization"""
    vec_sum_flops = algorithms["VecSum"].get_flops(input_length)
    vec_sum_err_branch_flops = algorithms["VecSumErrBranch"].get_flops(input_length)
    vec_sum_err_flops = algorithms["VecSumErr"].get_flops(
        (output_length * (output_length + 1) / 2.0)
    )

    return vec_sum_flops + vec_sum_err_branch_flops + vec_sum_err_flops


def get_multiplication_flops(expansion_length: int):
    """Calculate the number of floating point operations done by Multiplication"""
    return (
        expansion_length * (expansion_length + 1) / 2 * TWO_MULT_FMA_FLOPS
        + (expansion_length / 2) * algorithms["VecSum"].get_flops(expansion_length)
        + (expansion_length - 1) * 2
        + pow(expansion_length, 2)
        + get_renormalization_flops(expansion_length + 1, expansion_length)
    )


def get_fast_renormalization_flops(input_length: int, output_length: int):
    """Calculate the number of flops done by the reference Multiplication"""
    vec_sum_flops = algorithms["VecSumReference"].get_flops(input_length)
    vec_sum_err_branch_flops = (input_length - 2) * FAST_TWO_SUM_FLOPS
    vec_sum_err_flops = algorithms["VecSumErrReference"].get_flops(
        (output_length * (output_length + 1) / 2.0)
    )

    return vec_sum_flops + vec_sum_err_branch_flops + vec_sum_err_flops


def get_addition_reference_flops(a_length: int, b_length: int):
    """Calculate the number of flops done by the reference Addition"""
    expansion_length = a_length + b_length

    vec_sum_flops = algorithms["VecSum"].get_flops(expansion_length)
    vec_sum_err_branch_flops = (expansion_length - 2) * FAST_TWO_SUM_FLOPS

    return vec_sum_flops + vec_sum_err_branch_flops


# def get_multiplication_reference_flops(expansion_length: int):
#     """Calculate the number of flops done by the reference Multiplication"""
#     return (
#         expansion_length * (expansion_length + 1) / 2 * TWO_MULT_FMA_FLOPS
#         + (expansion_length / 2) * algorithms["VecSum"].get_flops(expansion_length)
#         + (expansion_length - 1) * 2
#         + pow(expansion_length, 2)
#         + get_renormalization_flops(expansion_length + 1, expansion_length)
#     )


def get_multiplication2_flops(expansion_length: int):
                              
    frexp_flops = 1
    ldexp_flops = 1

    mult_flops = 1
    add_flops = 1
    deposit_flops = TWO_SUM_FLOPS +1
    accumulate_flops = 2*deposit_flops 
    len_B = floor(expansion_length * 53 /45+2)
    return 2* len_B* frexp_flops + 2* ldexp_flops +(len_B-1)*mult_flops+ pow(expansion_length, 2) / 2 * (TWO_MULT_FMA_FLOPS  + mult_flops + accumulate_flops) + (expansion_length-1)*(2*mult_flops +  2*deposit_flops)+len_B*add_flops+ expansion_length*TWO_SUM_FLOPS

