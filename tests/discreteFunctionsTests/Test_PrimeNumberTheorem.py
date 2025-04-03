import unittest
from hypothesis import given, assume, settings, strategies as st
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))

from DiscreteFunctions import PrimeNumberTheorem

#python -m unittest tests/discreteFunctionsTests/Test_PrimeNumberTheorem.py

class TestPrimeNumberTheorem(unittest.TestCase):
    @given(st.integers(min_value=1000, max_value=10000))
    def test_pi_monotonicity(self, x):
        pi_x = PrimeNumberTheorem.π(x)
        pi_x_plus = PrimeNumberTheorem.π(x + 200)
        self.assertLessEqual(pi_x, pi_x_plus,
                             msg=f"π({x})={pi_x} > π({x+50})={pi_x_plus}")

    @given(st.integers(min_value=1000, max_value=10000))
    def test_continuous_estimate_bounds(self, x):
        pi_x = PrimeNumberTheorem.π(x)
        approx = PrimeNumberTheorem.approx_pnt_prime_count(x)
        self.assertAlmostEqual(approx, pi_x, delta=0.3 * pi_x,
                               msg=f"Continuous PNT estimate for {x} deviates: π(x)={pi_x}, approx={approx}")

    @given(st.integers(min_value=1000, max_value=10000))
    def test_summation_estimate_monotonic(self, x):
        est1 = PrimeNumberTheorem.pnt_prime_count(x)
        est2 = PrimeNumberTheorem.pnt_prime_count(x + 50)
        self.assertLessEqual(est1, est2,
                             msg=f"Discrete PNT estimate not monotonic: {est1} vs {est2}")

    @given(st.integers(min_value=1000, max_value=10000))
    def test_probability_bounds(self, n):
        p = PrimeNumberTheorem.prob(n)
        self.assertGreaterEqual(p, 0.0, msg=f"Probability for {n} below 0")
        self.assertLessEqual(p, 1.0, msg=f"Probability for {n} above 1")

if __name__ == '__main__':
    unittest.main()
