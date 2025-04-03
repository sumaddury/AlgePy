import unittest
from hypothesis import given, settings, strategies as st, assume
from fractions import Fraction
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))

from QuadraticStructures import Q
from SingletonStructures import Z

#python -m unittest tests/quadraticStructuresTests/Test_Q.py

class TestQ(unittest.TestCase):
    @given(st.integers(min_value=-1000, max_value=1000),
           st.integers(min_value=1, max_value=1000))
    def test_reduction(self, n, d):
        q = Q(n, d)
        frac = Fraction(n, d)
        self.assertEqual(int(q.n), frac.numerator,
                         msg=f"Numerator reduction error for {n}/{d}")
        self.assertEqual(int(q.d), frac.denominator,
                         msg=f"Denominator reduction error for {n}/{d}")

    @given(st.integers(min_value=-1000, max_value=1000),
           st.integers(min_value=1, max_value=1000))
    def test_addition_commutativity(self, n, d):
        q1 = Q(n, d)
        q2 = Q(n + 7, d)
        self.assertEqual(q1 + q2, q2 + q1,
                         msg="Addition in Q is not commutative")

    @given(st.integers(min_value=-1000, max_value=1000),
           st.integers(min_value=1, max_value=1000))
    def test_multiplication_commutativity(self, n, d):
        q1 = Q(n, d)
        q2 = Q(n + 3, d)
        self.assertEqual(q1 * q2, q2 * q1,
                         msg="Multiplication in Q is not commutative")

    @given(st.integers(min_value=-1000, max_value=1000),
           st.integers(min_value=1, max_value=1000))
    def test_reciprocal(self, n, d):
        assume(n != 0)
        q = Q(n, d)
        inv = q.reciprocal()
        prod = q * inv
        self.assertEqual(prod, Q.mul_id(),
                         msg=f"Reciprocal failed for {q}")

    def test_zero_and_mul_id(self):
        self.assertEqual(Q.zero(), Q(0, 1),
                         msg="Q.zero() not equal to 0/1")
        self.assertEqual(Q.mul_id(), Q(1, 1),
                         msg="Q.mul_id() not equal to 1/1")

if __name__ == '__main__':
    unittest.main()
