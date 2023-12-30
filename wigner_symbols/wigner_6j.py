import math
from typing import Tuple

from wigner_symbols.math_utils import factorial
from wigner_symbols.math_utils import get_triangle_coefficient
from wigner_symbols.math_utils import satisfies_trianglar_inequality


def calculate_6j(
    j1: float,
    j2: float,
    j3: float,
    j4: float,
    j5: float,
    j6: float,
) -> float:
    """Calculate the wigner 6j symbol.

    The 6j symbol is given as:
        { j1 j2 j3 }
        { j4 j5 j6 }
    where the explicit expression is given by:
        https://mathworld.wolfram.com/Wigner6j-Symbol.html

    args:
        j1
        j2
        j3
        j4
        j5
        j6

    returns:
        float: the wigner 6j symbol

    raise:
        NotAnIntegerException: 
            if any combination of j1, j2, j3, j4, j5, and j6,
            when summed/minused do not return an integer.
    """
    if not is_non_zero_symbol(j1, j2, j3, j4, j5, j6):
        return 0.0

    return math.prod(
        [
            get_square_root_term(j1, j2, j3, j4, j5, j6),
            get_summation_term(j1, j2, j3, j4, j5, j6),
        ],
    )


def get_square_root_term(
    j1: float,
    j2: float,
    j3: float,
    j4: float,
    j5: float,
    j6: float,
) -> float:
    return math.sqrt(
        math.prod(
            [
                get_triangle_coefficient(j1, j2, j3),
                get_triangle_coefficient(j1, j5, j6),
                get_triangle_coefficient(j4, j2, j6),
                get_triangle_coefficient(j4, j5, j3),
            ]
        )
    )


def get_summation_term(
    j1: float,
    j2: float,
    j3: float,
    j4: float,
    j5: float,
    j6: float,
) -> float:
    value = 0.0

    lower_bound, upper_bound = get_summation_limits(j1, j2, j3, j4, j5, j6)
    for k in range(lower_bound, upper_bound + 1):
        denominator = math.prod(
            [
                factorial(k - j1 - j2 - j3),
                factorial(k - j1 - j5 - j6),
                factorial(k - j4 - j2 - j6),
                factorial(k - j4 - j5 - j3),
                factorial(j1 + j2 + j4 + j5 - k),
                factorial(j2 + j3 + j5 + j6 - k),
                factorial(j3 + j1 + j6 + j4 - k),
            ],
        )
        value += ((-1.0) ** k) * factorial(k + 1) / denominator
    return value


def get_summation_limits(
    j1: float,
    j2: float,
    j3: float,
    j4: float,
    j5: float,
    j6: float,
) -> Tuple[int, int]:
    """Get the summation lower and upper bounds.

    Reference:
        1H. T. Johansson and C. Forssén,
        SIAM Journal on Scientific Computing, 2016, 38.
    
    Returns
        tuple[int, int]: the lower and upper bound values
            for carrying out the summation term in the
            calculation of the Wigner 6j symbol.
    """
    return (
        int(max(j1 + j2 + j3, j4 + j5 + j3, j1 + j5 + j6, j4 + j2 + j6)),
        int(min(j1 + j2 + j4 + j5, j1 + j3 + j4 + j6, j2 + j3 + j5 + j6)),
    )


def is_non_zero_symbol(
    j1: float, j2: float, j3: float, j4: float, j5: float, j6: float
) -> bool:
    if (
        satisfies_trianglar_inequality(j1, j2, j3)
        and satisfies_trianglar_inequality(j1, j5, j6)
        and satisfies_trianglar_inequality(j4, j2, j6)
        and satisfies_trianglar_inequality(j4, j5, j3)
    ):
        return True
    else:
        return False
