import unittest
from hypothesis import given, assume, settings, strategies as st
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from .QuadraticStructures import QuadRatField, QuadRat, Q, QuadIntRing
from .SingletonStructures import Z
from .DiscreteFunctions import ArithmeticFunctions

#python -m unittest tests/quadraticStructuresTests/Test_QuadIntRing.py

class TestQuadRatField(unittest.TestCase):
    d = 2

    def setUp(self):
        self.Field = QuadRatField(self.d)

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_field_additive_inverse(self, n1, d1, n2, d2):
        x = self.Field(Q(n1, d1), Q(n2, d2))
        zero = self.Field(Q(0, 1), Q(0, 1))
        self.assertEqual(x + (-x), zero, "x + (-x) must equal zero in the field")

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_field_multiplicative_inverse(self, n1, d1, n2, d2):
        x = self.Field(Q(n1, d1), Q(n2, d2))
        zero = self.Field(Q(0, 1), Q(0, 1))
        assume(x != zero)
        one = self.Field(Q.mul_id(), Q.zero())
        self.assertEqual(x * x.inverse(), one, "x * x.inverse() must equal one in the field")

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_field_distributivity(self, n1, d1, n2, d2, n3, d3, n4, d4):
        x = self.Field(Q(n1, d1), Q(n2, d2))
        y = self.Field(Q(n3, d3), Q(n4, d4))
        u = self.Field(Q(random.randint(-100, 100), 1), Q(random.randint(-100, 100), 1))
        self.assertEqual(u * (x + y), u * x + u * y, "Distributivity must hold in the field")

    def test_integer_ring(self):
        IR1 = self.Field.integer_ring()
        IR2 = QuadIntRing(self.d, force_ufd=True)
        for a in range(-50, 51, 10):
            for b in range(-50, 51, 10):
                self.assertEqual(IR1(a, b), IR2(a, b), "Integer ring construction must be consistent")

    def test_discriminant(self):
        disc = self.Field.discriminant()
        expected = Z(4) * Z(self.d)
        self.assertEqual(disc, expected, "Field discriminant is incorrect for d=2")

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_embeddings(self, n1, d1, n2, d2):
        x = self.Field(Q(n1, d1), Q(n2, d2))
        emb = self.Field.embed(x)
        emb_prime = self.Field.embed_prime(x)
        rpart = float(x.a)
        self.assertAlmostEqual(emb + emb_prime, 2 * rpart, delta=0.001,
                               msg="Sum of embeddings must equal twice the rational part")

    @given(
        st.integers(min_value=-10, max_value=10)
    )
    def test_exponentiation(self, exp):
        a_num = random.randint(-1000, 1000)
        a_den = random.randint(1, 1000)
        b_num = random.randint(-1000, 1000)
        b_den = random.randint(1, 1000)
        x = self.Field(Q(a_num, a_den), Q(b_num, b_den))
        power = x ** Z(exp)
        prod = self.Field(Q.mul_id(), Q.zero())
        for _ in range(exp):
            prod = prod * x
        self.assertEqual(power, prod, "Exponentiation must equal repeated multiplication")

    @given(
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000),
        st.integers(min_value=-1000, max_value=1000), st.integers(min_value=1, max_value=1000)
    )
    def test_inversion_consistency(self, n1, d1, n2, d2):
        x = self.Field(Q(n1, d1), Q(n2, d2))
        zero = self.Field(Q(0, 1), Q(0, 1))
        assume(x != zero)
        one = self.Field(Q.mul_id(), Q.zero())
        self.assertEqual(x * x.inverse(), one, "Inversion must be consistent in the field")

if __name__ == '__main__':
    unittest.main()
