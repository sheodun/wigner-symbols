import math
import time

import pytest

from wigner_symbols import constants
from wigner_symbols import calculate_3j


class TestWigner3j:

    @pytest.mark.parametrize(
        "jm, expected_value",
        [
            ((2.0, 2.0, 2.0, 0.0, 0.0, 0.0), -math.sqrt(2 / 35)),
            ((2.0, 2.0, 2.0, 2.0, 0.0, -2.0), math.sqrt(2 / 35)),
            ((6.0, 4.0, 2.0, 0.0, 0.0, 0.0), math.sqrt(5 / 143)),
            ((6.0, 4.0, 2.0, 1.0, 0.0, 0.0), 0.0),
            ((2.5, 1.5, 1.0, 1.5, -0.5, -1.0), -math.sqrt(1 / 10)),
        ],
    )
    def test_3j_calculation(self, jm, expected_value):
        start = time.time()
        output = calculate_3j(*jm)
        finish = time.time()

        assert abs(expected_value - output) < constants.EPSILON
