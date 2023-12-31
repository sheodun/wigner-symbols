import math
import time

import pytest

from wigner_symbols import constants
from wigner_symbols import calculate_6j


class TestWigner6j:

    @pytest.mark.parametrize(
        "jm, expected_value",
        [
            ((2.0, 2.0, 2.0, 2.0, 2.0, 2.0), -3.0 / 70.0),
            ((1.0, 2.0, 3.0, 2.0, 1.0, 2.0), 1.0 / (5.0 * math.sqrt(21.0))),
            ((0.5, 1.5, 4.0, 3.0, 2.0, 0.5), 0.0),
        ],
    )
    def test_6j_calculation(self, jm, expected_value):
        start = time.time()
        output = calculate_6j(*jm)
        finish = time.time()

        assert abs(expected_value - output) < constants.EPSILON
