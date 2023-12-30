import math
from typing import Tuple

from wigner_symbols import constants
from wigner_symbols.math_utils import factorial
from wigner_symbols.math_utils import get_triangle_coefficient
from wigner_symbols.math_utils import satisfies_trianglar_inequality


def calculate_3j(
    j1: float,
    j2: float,
    j3: float,
    m1: float,
    m2: float,
    m3: float,
) -> float:
    """Calculate the wigner 3j symbol.

    The 3j symbol is given as:
        | j1 j2 j3 |
        | m1 m2 m3 |
    where the explicit calculation is given:
        https://en.wikipedia.org/wiki/3-j_symbol#Explicit_expression

    args:
        j1
        j2
        j3
        m1
        m2
        m3

    returns:
        float: the wigner 3j symbol

    raise:
        NotAnIntegerException: 
            if any combination of j1, j2, j3, m1, m2, and m3,
            when summed/minused do not return an integer.
    """
    if not is_non_zero_symbol(j1, j2, j3, m1, m2, m3):
        return 0.0

    phase = get_phase(j1, j2, m3)
    return math.prod(
        [
            phase,
            math.sqrt(get_triangle_coefficient(j1, j2, j3)),
            get_square_root_factorial_term(j1, j2, j3, m1, m2, m3),
            get_summation_term(j1, j2, j3, m1, m2),
        ]
    )


def get_summation_term(
    j1: float,
    j2: float,
    j3: float,
    m1: float,
    m2: float,
) -> float:
    """Return the summation term for calculating the wigner 3j symbol."""
    value = 0.0
    lower_bound, upper_bound = get_summation_limits(j1, j2, j3, m1, m2)
    for k in range(lower_bound, upper_bound + 1):
        denominator = math.prod(
            [
                factorial(k),
                factorial(j1 + j2 - j3 - k),
                factorial(j1 - m1 - k),
                factorial(j2 + m2 - k),
                factorial(j3 - j2 + m1 + k),
                factorial(j3 - j1 - m2 + k),
            ],
        )
        value += ((-1.0) ** k) / denominator
    return value


def get_summation_limits(
    j1: float,
    j2: float,
    j3: float,
    m1: float,
    m2: float,
) -> Tuple[int, int]:
    """Get the summation lower and upper bounds.
    
    Returns
        tuple[int, int]: the lower and upper bound values
            for carrying out the summation term in the
            calculation of the Wigner 3j symbol.
    """
    return (
        max(0, max(int(j2 - j3 - m1), int(j1 - j3 + m2))),
        min(
            int(j1 + j2 - j3),
            min(
                int(j1 - m1),
                int(j2 + m2),
            ),
        ),
    )


def get_square_root_factorial_term(
    j1: float, j2: float, j3: float, m1: float, m2: float, m3: float
) -> float:
    """Return the square root factorial term for calculating 3j symbol."""
    return math.sqrt(
        math.prod(
            [
                factorial(j1 - m1),
                factorial(j1 + m1),
                factorial(j2 - m2),
                factorial(j2 + m2),
                factorial(j3 - m3),
                factorial(j3 + m3),
            ]
        )
    )


def get_phase(j1: float, j2: float, m3: float) -> float:
    return (-1.0) ** (j1 - j2 - m3)


def is_non_zero_symbol(
    j1: float, j2: float, j3: float, m1: float, m2: float, m3: float
) -> bool:
    """
    Check if the combination of angular momentum numbers means
    that the 3j symbol is equal to zero by defined selection rules.
    """
    if not kronecker_delta(m1, m2, m3):
        return False

    for j, m in ((j1, m1), (j2, m2), (j3, m3)):
        if not is_projection_within_range(j, m):
            return False

    if not satisfies_trianglar_inequality(j1, j2, j3):
        return False

    if not is_sum_of_j_an_integer(j1, j2, j3):
        return False

    return True


def is_sum_of_j_an_integer(j1: float, j2: float, j3: float) -> bool:
    return abs(j1 + j2 + j3 - int(j1 + j2 + j3)) < constants.EPSILON


def is_projection_within_range(j: float, m: float) -> bool:
    return -j <= m <= j


def kronecker_delta(m1: float, m2: float, m3: float) -> bool:
    return abs(m1 + m2 + m3) < constants.EPSILON
