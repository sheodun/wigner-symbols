import math

from wigner_symbols import constants
from wigner_symbols.exceptions import NotAnIntegerException


def factorial(value: float) -> float:
    """Utility method to calculate factorial.

    Args:
        value: the value to calculate the factorial for.

    Raises:
        NotAnIntegerException: if the float passed through
            is not close to the integer value
    """
    integer_value = int(value)
    if abs(value - integer_value) > constants.EPSILON:
        raise NotAnIntegerException(
            f"the float number {value} is not a floating integer."
        )
    return math.factorial(integer_value)


def get_triangle_coefficient(j1: float, j2: float, j3: float) -> float:
    """Calculate the triangle coefficient for three angular momenta.

    https://mathworld.wolfram.com/TriangleCoefficient.html

    Args:
        j1: the first angular momentum number. Either integer or
            half-integer.
        j2: the second angular momentum number. Either integer or
            half-integer.
        j3: the third angular momentum number. Either integer or
            half-integer.

    Returns:
        float: the triangle coefficient.
    """
    numerator = math.prod(
        [
            factorial(j1 + j2 - j3),
            factorial(j1 - j2 + j3),
            factorial(-j1 + j2 + j3),
        ]
    )
    return numerator / factorial(j1 + j2 + j3 + 1)


def satisfies_trianglar_inequality(j1: float, j2: float, j3: float) -> bool:
    return abs(j1 - j2) <= j3 <= j1 + j2
