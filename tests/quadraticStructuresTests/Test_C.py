import unittest
from hypothesis import given, settings, strategies as st
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from QuadraticStructures import C
from SingletonStructures import R, Z

#python -m unittest tests/quadraticStructuresTests/Test_C.py

class TestC(unittest.TestCase):
    @given(st.floats(min_value=-1000, max_value=1000, allow_nan=False, allow_infinity=False),
           st.floats(min_value=-1000, max_value=1000, allow_nan=False, allow_infinity=False))
    def test_addition_commutativity(self, a, b):
        c1 = C(a, b)
        c2 = C(b, a)
        res1 = c1 + c2
        res2 = c2 + c1
        self.assertAlmostEqual(res1.to_tuple()[0], res2.to_tuple()[0], places=6,
                               msg=f"Addition (real part) not commutative for {c1} and {c2}")
        self.assertAlmostEqual(res1.to_tuple()[1], res2.to_tuple()[1], places=6,
                               msg=f"Addition (imaginary part) not commutative for {c1} and {c2}")

    @given(st.floats(min_value=-500, max_value=500, allow_nan=False, allow_infinity=False),
           st.floats(min_value=-500, max_value=500, allow_nan=False, allow_infinity=False))
    def test_multiplication_associativity(self, a, b):
        c1 = C(a, b)
        c2 = C(b, a)
        c3 = C(a + 1, b + 1)
        res1 = (c1 * c2) * c3
        res2 = c1 * (c2 * c3)
        self.assertAlmostEqual(res1.to_tuple()[0], res2.to_tuple()[0], places=6,
                               msg="Multiplication (real part) not associative")
        self.assertAlmostEqual(res1.to_tuple()[1], res2.to_tuple()[1], places=6,
                               msg="Multiplication (imaginary part) not associative")

    @given(st.floats(min_value=-1000, max_value=1000, allow_nan=False, allow_infinity=False),
           st.floats(min_value=-1000, max_value=1000, allow_nan=False, allow_infinity=False))
    def test_conjugate_norm(self, a, b):
        c_obj = C(a, b)
        prod = c_obj * c_obj.conjugate()
        self.assertAlmostEqual(c_obj.conjugate().to_tuple()[1], -c_obj.to_tuple()[1], places=6,
                               msg="Conjugation did not negate the imaginary part")
        expected_norm = math.sqrt(a**2 + b**2)
        self.assertAlmostEqual(float(C.norm(c_obj).r), expected_norm, places=6,
                               msg=f"Norm mismatch for {c_obj}")

    @given(st.integers(min_value=0, max_value=10))
    def test_power(self, exp):
        base = C(3.0, 4.0)
        result = base ** Z(exp)
        prod = C.mul_id()
        for _ in range(exp):
            prod = prod * base
        self.assertAlmostEqual(result.to_tuple()[0], prod.to_tuple()[0], places=6,
                               msg=f"Power real part failed for exponent {exp}")
        self.assertAlmostEqual(result.to_tuple()[1], prod.to_tuple()[1], places=6,
                               msg=f"Power imaginary part failed for exponent {exp}")

    def test_zero_and_mul_id(self):
        self.assertEqual(C.zero(), C(0, 0),
                         msg="C.zero() does not return 0+0i")
        self.assertEqual(C.mul_id(), C(1, 0),
                         msg="C.mul_id() does not return 1+0i")

if __name__ == '__main__':
    unittest.main()