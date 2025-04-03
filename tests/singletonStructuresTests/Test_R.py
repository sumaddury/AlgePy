import unittest
from hypothesis import given, assume, settings, strategies as st
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from .SingletonStructures import R

#python -m unittest tests/singletonStructuresTests/Test_R.py

class TestR(unittest.TestCase):

    @given(st.floats(min_value=-1000, max_value=1000, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-1000, max_value=1000, allow_infinity=False, allow_nan=False))
    def test_add_commutative(self, a, b):
        A, B = R(a), R(b)
        left = A + B
        right = B + A
        self.assertAlmostEqual(float(left), float(right), places=7,
                               msg=f"Addition not commutative for {a} + {b}")

    @given(st.floats(min_value=-100, max_value=100, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-100, max_value=100, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-100, max_value=100, allow_infinity=False, allow_nan=False))
    def test_add_associative(self, a, b, c):
        A, B, C = R(a), R(b), R(c)
        left = (A + B) + C
        right = A + (B + C)
        self.assertAlmostEqual(float(left), float(right), places=7,
                               msg=f"Addition not associative for {a},{b},{c}")

    @given(st.floats(min_value=-100, max_value=100, allow_infinity=False, allow_nan=False))
    def test_additive_identity(self, a):
        A = R(a)
        zero = R.zero()
        self.assertAlmostEqual(float(A + zero), a, places=7,
                               msg=f"Additive identity failed for {a}")

    @given(st.floats(min_value=-100, max_value=100, allow_infinity=False, allow_nan=False))
    def test_negation(self, a):
        A = R(a)
        negA = -A
        res = A + negA
        self.assertAlmostEqual(float(res), 0.0, places=7,
                               msg=f"Additive inverse failed for {a}")

    @given(st.floats(min_value=-50, max_value=50, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-50, max_value=50, allow_infinity=False, allow_nan=False))
    def test_mul_commutative(self, a, b):
        A, B = R(a), R(b)
        left = A * B
        right = B * A
        self.assertAlmostEqual(float(left), float(right), places=7,
                               msg=f"Multiplication not commutative for {a} * {b}")

    @given(st.floats(min_value=-20, max_value=20, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-20, max_value=20, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-20, max_value=20, allow_infinity=False, allow_nan=False))
    def test_mul_associative(self, a, b, c):
        A, B, C = R(a), R(b), R(c)
        left = (A * B) * C
        right = A * (B * C)
        self.assertAlmostEqual(float(left), float(right), places=7,
                               msg=f"Multiplication not associative for {a},{b},{c}")

    @given(st.floats(min_value=-20, max_value=20, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-20, max_value=20, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-20, max_value=20, allow_infinity=False, allow_nan=False))
    def test_distributive(self, a, b, c):
        A, B, C = R(a), R(b), R(c)
        left = A * (B + C)
        right = (A * B) + (A * C)
        self.assertAlmostEqual(float(left), float(right), places=7,
                               msg=f"Distributivity failed for {a},{b},{c}")

    @given(st.floats(min_value=-20, max_value=20, allow_infinity=False, allow_nan=False))
    def test_mul_identity(self, a):
        A = R(a)
        one = R.mul_id()
        res = A * one
        self.assertAlmostEqual(float(res), a, places=7,
                               msg=f"Multiplicative identity failed for {a}")

    @given(st.floats(min_value=1e-5, max_value=20, allow_infinity=False, allow_nan=False))
    def test_division_inverse(self, a):
        """
        For a != 0, a / a == 1.0
        """
        A = R(a)
        inv = A / A
        self.assertAlmostEqual(float(inv), 1.0, places=7,
                               msg=f"Division inverse failed for {a}")

    @given(st.floats(min_value=-5, max_value=5, allow_infinity=False, allow_nan=False),
           st.floats(min_value=-2, max_value=2, allow_infinity=False, allow_nan=False))
    def test_pow(self, base, exponent):
        """
        Check exponentiation in R with float exponent. 
        We skip negative base with non-integer exponent to avoid domain errors.
        """
        assume(base != 0)
        if base < 0 and not exponent.is_integer():
            return
        A = R(base)
        E = R(exponent)
        try:
            result = A ** E
        except ValueError:
            return
        expected = float(base) ** float(exponent)
        self.assertAlmostEqual(float(result), expected, places=5,
                               msg=f"Power mismatch for {base}^{exponent}")

    @given(st.floats(min_value=-100, max_value=100, allow_infinity=False, allow_nan=False))
    def test_abs(self, a):
        A = R(a)
        absA = R.abs(A)
        self.assertGreaterEqual(float(absA), 0.0,
                                msg=f"Absolute value < 0 ??? for {a}")
        self.assertAlmostEqual(float(absA), abs(a), places=7,
                               msg=f"abs mismatch for {a}")

if __name__ == "__main__":
    unittest.main()
