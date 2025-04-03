import unittest
from hypothesis import given, assume, settings, strategies as st
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from .QuadraticStructures import QuadRat, Q
from .SingletonStructures import Z

#python -m unittest tests/quadraticStructuresTests/Test_QuadRat.py

class TestQuadRat(unittest.TestCase):
    d = 2

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_addition_commutativity(self, n1, d1, n2, d2):
        a1 = Q(n1, d1)
        a2 = Q(n2, d2)
        x = QuadRat(a1, a2, self.d)
        
        y = QuadRat(a2, a1, self.d)
        
        self.assertEqual(x + y, y + x, "Addition must be commutative")

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_subtraction_and_negation(self, n1, d1, n2, d2):
        a1 = Q(n1, d1)
        a2 = Q(n2, d2)
        x = QuadRat(a1, a2, self.d)
        zero = QuadRat(Q(0,1), Q(0,1), self.d)
        self.assertEqual(x - x, zero, "x - x must equal zero")
        self.assertEqual(x + (-x), zero, "x + (-x) must equal zero")

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_multiplicative_inverse(self, n1, d1, n2, d2):
        a1 = Q(n1, d1)
        a2 = Q(n2, d2)
        x = QuadRat(a1, a2, self.d)
        assume(not QuadRat.norm(x) == Q.zero())
        zero = QuadRat(Q(0,1), Q(0,1), self.d)
        assume(x != zero)
        one = QuadRat(Q.mul_id(), Q.zero(), self.d)
        self.assertEqual(x * x.inverse(), one, "x * x.inverse() must equal one")

    @given(
        st.integers(min_value=-100, max_value=100), st.integers(min_value=1, max_value=100),
        st.integers(min_value=-100, max_value=100), st.integers(min_value=1, max_value=100),
        st.integers(min_value=-100, max_value=100), st.integers(min_value=1, max_value=100),
        st.integers(min_value=-100, max_value=100), st.integers(min_value=1, max_value=100)
    )
    def test_distributivity(self, n1, d1, n2, d2, n3, d3, n4, d4):
        a1 = Q(n1, d1)
        a2 = Q(n2, d2)
        b1 = Q(n3, d3)
        b2 = Q(n4, d4)
        x = QuadRat(a1, a2, self.d)
        y = QuadRat(b1, b2, self.d)
        u = QuadRat(Q(random.randint(-100,100), 1), Q(random.randint(-100,100), 1), self.d)
        self.assertEqual(u * (x + y), u*x + u*y, "Distributivity must hold")

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_trace(self, n1, d1, n2, d2):
        a1 = Q(n1, d1)
        a2 = Q(n2, d2)
        x = QuadRat(a1, a2, self.d)
        expected = (x + x.conjugate()).a
        self.assertEqual(QuadRat.trace(x), expected, "Trace must equal x + conjugate(x) (rational part)")

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_conjugate_norm_relationship(self, n1, d1, n2, d2):
        a1 = Q(n1, d1)
        a2 = Q(n2, d2)
        x = QuadRat(a1, a2, self.d)
        prod = x * x.conjugate()
        norm_val = QuadRat.norm(x)
        expected = QuadRat(norm_val, Q.zero(), self.d)
        self.assertEqual(prod, expected, "x * conjugate(x) must equal norm(x) as a field element")

    @given(
        st.integers(min_value=0, max_value=100)
    )
    def test_exponentiation(self, exp):
        a_num = random.randint(-100, 100)
        a_den = random.randint(1, 100)
        b_num = random.randint(-100, 100)
        b_den = random.randint(1, 100)
        a = Q(a_num, a_den)
        b = Q(b_num, b_den)
        x = QuadRat(a, b, self.d)
        power = x ** Z(exp)
        prod = QuadRat(Q.mul_id(), Q.zero(), self.d)
        for _ in range(exp):
            prod = prod * x
        self.assertEqual(power, prod, "Exponentiation must match repeated multiplication")

if __name__ == '__main__':
    unittest.main()
