import unittest
from hypothesis import given, assume, settings, strategies as st
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from .QuadraticStructures import QuadInt
from .SingletonStructures import Z

#python -m unittest tests/quadraticStructuresTests/Test_QuadInt.py

class TestQuadIntStandard(unittest.TestCase):
    def setUp(self):
        self.d = 2
    
    @given(st.integers(min_value=-50, max_value=50),
           st.integers(min_value=-50, max_value=50))
    def test_norm_standard(self, a, b):
        x = QuadInt(Z(a), Z(b), Z(self.d))
        product = x * x.conjugate()
        self.assertEqual(int(product.b), 0, msg=f"Standard: Secondary part nonzero for {x}")
        expected = a * a - self.d * (b * b)
        self.assertEqual(int(product.a), expected, msg=f"Standard: Norm mismatch for {x}: got {int(product.a)}, expected {expected}")
    
    @given(st.integers(min_value=-50, max_value=50),
           st.integers(min_value=-50, max_value=50))
    def test_addition_multiplication(self, a, b):
        x = QuadInt(Z(a), Z(b), Z(self.d))
        y = QuadInt(Z(b), Z(a), Z(self.d))
        self.assertEqual(x + y, y + x, msg="Standard: Addition not commutative")
        self.assertEqual(x * y, y * x, msg="Standard: Multiplication not commutative")
    
    @given(st.integers(min_value=-10, max_value=10))
    def test_inverse_unit_standard(self, a):
        assume(abs(a) == 1)
        x = QuadInt(Z(a), Z(0), Z(self.d))
        inv = x.inverse()
        self.assertEqual(x * inv, QuadInt(1,0,self.d), msg=f"Standard: Inverse failed for unit {x}")
    
    @given(st.integers(min_value=-20, max_value=20),
           st.integers(min_value=-20, max_value=20),
           st.integers(min_value=0, max_value=10))
    def test_power_standard(self, a, b, exp):
        x = QuadInt(Z(a), Z(b), Z(self.d))
        prod = QuadInt(1,0,self.d)
        for _ in range(exp):
            prod = prod * x
        self.assertEqual(x ** Z(exp), prod, msg=f"Standard: Power function failed for {x} with exponent {exp}")
    
    def test_identity_standard(self):
        self.assertEqual(QuadInt(0,0,self.d), QuadInt(Z(0), Z(0), Z(self.d)),
                         msg="Standard: zero() incorrect")
        self.assertEqual(QuadInt(1,0,self.d), QuadInt(Z(1), Z(0), Z(self.d)),
                         msg="Standard: mul_id() incorrect")

class TestQuadIntHalfIntegral(unittest.TestCase):
    """
    Tests for QuadInt using the half-integral representation.
    Here we use d = 5 (since 5 ≡ 1 mod 4) so that elements are expressed canonically
    as a + bω with ω = (1+√5)/2 and the norm is: N(x)=a² + ab - b².
    """
    def setUp(self):
        self.d = 5
    
    @given(st.integers(min_value=-50, max_value=50),
           st.integers(min_value=-50, max_value=50))
    def test_norm_half_integral(self, a, b):
        x = QuadInt(Z(a), Z(b), Z(self.d))
        product = x * x.conjugate()
        self.assertEqual(int(product.b), 0, msg=f"Half-integral: Secondary part nonzero for {x}")
        expected = a * a + a * b - b * b
        self.assertEqual(int(product.a), expected, msg=f"Half-integral: Norm mismatch for {x}: got {int(product.a)}, expected {expected}")
    
    @given(st.integers(min_value=-50, max_value=50),
           st.integers(min_value=-50, max_value=50))
    def test_addition_multiplication_half_integral(self, a, b):
        x = QuadInt(Z(a), Z(b), Z(self.d))
        y = QuadInt(Z(b), Z(a), Z(self.d))
        self.assertEqual(x + y, y + x, msg="Half-integral: Addition not commutative")
        self.assertEqual(x * y, y * x, msg="Half-integral: Multiplication not commutative")
    
    @given(st.integers(min_value=-10, max_value=10))
    def test_inverse_unit_half_integral(self, a):
        assume(abs(a) == 1)
        x = QuadInt(Z(a), Z(0), Z(self.d))
        inv = x.inverse()
        self.assertEqual(x * inv, QuadInt(1,0,self.d), msg=f"Half-integral: Inverse failed for unit {x}")
    
    @given(st.integers(min_value=-20, max_value=20),
           st.integers(min_value=-20, max_value=20),
           st.integers(min_value=0, max_value=10))
    def test_power_half_integral(self, a, b, exp):
        x = QuadInt(Z(a), Z(b), Z(self.d))
        prod = QuadInt(1,0,self.d)
        for _ in range(exp):
            prod = prod * x
        self.assertEqual(x ** Z(exp), prod, msg=f"Half-integral: Power function failed for {x} with exponent {exp}")
    
    def test_identity_half_integral(self):
        self.assertEqual(QuadInt(0,0,self.d), QuadInt(Z(0), Z(0), Z(self.d)),
                         msg="Half-integral: zero() incorrect")
        self.assertEqual(QuadInt(1,0,self.d), QuadInt(Z(1), Z(0), Z(self.d)),
                         msg="Half-integral: mul_id() incorrect")

if __name__ == '__main__':
    unittest.main()
