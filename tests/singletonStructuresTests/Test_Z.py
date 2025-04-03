import unittest
from hypothesis import given, assume, settings, strategies as st
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from .SingletonStructures import Z

#python -m unittest tests/singletonStructuresTests/Test_Z.py

class TestZ(unittest.TestCase):

    @given(st.integers(min_value=-1000, max_value=1000),
           st.integers(min_value=-1000, max_value=1000))
    def test_addition_commutative(self, a, b):
        """
        Commutativity: a + b = b + a
        """
        A, B = Z(a), Z(b)
        self.assertEqual(A + B, B + A,
                         msg=f"Addition not commutative for {a} and {b}")

    @given(st.integers(min_value=-500, max_value=500),
           st.integers(min_value=-500, max_value=500),
           st.integers(min_value=-500, max_value=500))
    def test_addition_associative(self, a, b, c):
        """
        Associativity: (a + b) + c = a + (b + c)
        """
        A, B, C = Z(a), Z(b), Z(c)
        left = (A + B) + C
        right = A + (B + C)
        self.assertEqual(left, right,
                         msg=f"Addition not associative for {a},{b},{c}")

    @given(st.integers(min_value=-500, max_value=500))
    def test_additive_identity(self, a):
        """
        Existence of additive identity: a + 0 = a
        """
        A = Z(a)
        self.assertEqual(A + Z.zero(), A,
                         msg=f"Additive identity failed for {a}")

    @given(st.integers(min_value=-500, max_value=500))
    def test_negation(self, a):
        """
        Existence of additive inverse: a + (-a) = 0
        """
        A = Z(a)
        self.assertEqual(A + (-A), Z.zero(),
                         msg=f"Additive inverse failed for {a}")

    @given(st.integers(min_value=-50, max_value=50),
           st.integers(min_value=-50, max_value=50))
    def test_multiplication_commutative(self, a, b):
        """
        Commutativity of multiplication: a * b = b * a
        """
        A, B = Z(a), Z(b)
        self.assertEqual(A * B, B * A,
                         msg=f"Multiplication not commutative for {a} and {b}")

    @given(st.integers(min_value=-20, max_value=20),
           st.integers(min_value=-20, max_value=20),
           st.integers(min_value=-20, max_value=20))
    def test_multiplication_associative(self, a, b, c):
        """
        Associativity of multiplication: (a*b)*c = a*(b*c)
        """
        A, B, C = Z(a), Z(b), Z(c)
        left = (A * B) * C
        right = A * (B * C)
        self.assertEqual(left, right,
                         msg=f"Multiplication not associative for {a},{b},{c}")

    @given(st.integers(min_value=-20, max_value=20),
           st.integers(min_value=-20, max_value=20),
           st.integers(min_value=-20, max_value=20))
    def test_distributive(self, a, b, c):
        """
        Distributivity: a*(b + c) = a*b + a*c
        """
        A, B, C = Z(a), Z(b), Z(c)
        left = A * (B + C)
        right = (A * B) + (A * C)
        self.assertEqual(left, right,
                         msg=f"Distributivity failed for {a},{b},{c}")

    @given(st.integers(min_value=-100, max_value=100))
    def test_multiplicative_identity(self, a):
        """
        Existence of multiplicative identity: a * 1 = a
        """
        A = Z(a)
        self.assertEqual(A * Z.mul_id(), A,
                         msg=f"Multiplicative identity failed for {a}")

    @given(st.integers(min_value=-20, max_value=20),
           st.integers(min_value=0, max_value=5))
    def test_power(self, a, exp):
        """
        Exponentiation: a^exp in Z
        """
        A = Z(a)
        E = Z(exp)
        result = A ** E
        expected = Z(a ** exp)
        self.assertEqual(result, expected,
                         msg=f"Power mismatch for {a}^{exp}")

    @given(st.integers(min_value=-50, max_value=50),
           st.integers(min_value=1, max_value=50))
    def test_floordiv_mod(self, a, b):
        """
        Floor division and mod: check (a // b)*b + (a % b) = a
        """
        assume(b != 0)
        A, B = Z(a), Z(b)
        q = A // B
        r = A % B
        lhs = (q * B) + r
        rhs = Z(a)
        self.assertEqual(lhs, rhs,
                         msg=f"Division algorithm property failed for {a}, {b}")

    @given(st.integers(min_value=-100, max_value=100))
    def test_abs(self, a):
        """
        Absolute value property: abs(a) >= 0 and is consistent with python's abs.
        """
        A = Z(a)
        absA = Z.abs(A)
        self.assertGreaterEqual(int(absA), 0,
                                msg=f"abs({a}) < 0 ???")
        self.assertEqual(int(absA), abs(a),
                         msg=f"abs mismatch for {a}")

if __name__ == "__main__":
    unittest.main()
