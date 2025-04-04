import unittest
from hypothesis import given, settings, strategies as st, assume
import math
import random
import sys
import os

# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from algepy.QuadraticStructures import QuadIntRing, QuadInt
from algepy.SingletonStructures import Z
from algepy.DiscreteFunctions import Factorization, ArithmeticFunctions

#python -m unittest tests/quadraticStructuresTests/Test_QuadIntRing.py

class TestQuadIntRingStandard(unittest.TestCase):
    def setUp(self):
        self.d = 2
        self.Ring = QuadIntRing(self.d, force_ufd=True)
    
    @given(st.integers(min_value=-100, max_value=100))
    def test_gcd_standard(self, a):
        x = self.Ring(a, 0)
        y = self.Ring(a + 14, 0)
        g = self.Ring.gcd(x, y)
        expected = math.gcd(a, a + 14)
        self.assertEqual(abs(int(g.a)), expected,
                         msg=f"Standard: gcd mismatch for {a} and {a+14}")
        
    def test_fundamental_unit_standard(self):
        """
        For d = 2, the canonical fundamental unit is 1 + √2 (up to unit multiplication).
        Its embedding should have absolute value approximately 1 + √2.
        """
        fu = self.Ring.get_fundamental_unit()
        emb = abs(float(self.Ring.embed(fu)))
        expected = 1 + math.sqrt(2)  # approximately 2.41421356...
        self.assertAlmostEqual(emb, expected, delta=0.001,
            msg=f"Standard: Fundamental unit embedding mismatch for d=2: got {emb}, expected {expected}")
        
    def test_imaginary_units_minus1(self):
        d = -1
        Ring = QuadIntRing(d, force_ufd=True)
        units = Ring.get_imaginary_units()
        # Compute embeddings as complex numbers
        embeddings = {complex(*Ring.embed(u).to_tuple()) for u in units}
        expected = {1+0j, -1+0j, 0+1j, 0-1j}
        self.assertEqual(embeddings, expected,
            msg=f"Imaginary: Units for d=-1 are incorrect: got {embeddings}")

        
    @given(st.integers(min_value=-50, max_value=50), st.integers(min_value=-50, max_value=50))
    def test_normalize_standard(self, a, b):
        assume(not (a == 0 and b == 0))
        z = self.Ring(a, b)
        zn = z.normalize()
        self.assertEqual(zn.normalize(), zn,
                         msg=f"Standard normalization not idempotent for element {z}")
        
    @given(st.integers(min_value=-50, max_value=50), st.integers(min_value=-50, max_value=50))
    def test_normalize_embedding_range_standard(self, a, b):
        assume(not (a == 0 and b == 0))
        z = self.Ring(a, b).normalize()
        if abs(int(QuadInt.norm(z))) == 1:
            return
        sigma_z = abs(float(self.Ring.embed(z)))
        epsilon = self.Ring.get_fundamental_unit()
        sigma_epsilon = abs(float(self.Ring.embed(epsilon)))
        self.assertGreaterEqual(sigma_z, 1,
                                msg=f"Standard: Embedding of normalized {z} is less than 1: {sigma_z}")
        self.assertLess(sigma_z, sigma_epsilon,
                        msg=f"Standard: Embedding of normalized {z} is not less than |σ(ε)| = {sigma_epsilon}")
    
    def test_normalize_idempotence_standard_imaginary(self):
        d = -2
        Ring = QuadIntRing(d, force_ufd=True)
        # Use a few fixed test cases.
        test_elements = [
            Ring(3, 4),    # 3 + 4√(-2)
            Ring(-5, 2),   # -5 + 2√(-2)
            Ring(7, -3)    # 7 - 3√(-2)
        ]
        for z in test_elements:
            zn = z.normalize()
            self.assertEqual(zn.normalize(), zn,
                             msg=f"Standard imaginary normalization not idempotent for element {z}")
            
    def test_normalize_embedding_minimal_standard_imaginary(self):
        d = -2
        Ring = QuadIntRing(d, force_ufd=True)
        z = Ring(3, 4)  # Example element: 3 + 4√(-2)
        normed = z.normalize()
        # The only associates are z and -z.
        candidate1 = z
        candidate2 = -z
        emb1 = Ring.embed(candidate1).to_tuple()
        emb2 = Ring.embed(candidate2).to_tuple()
        best = min(emb1, emb2)
        normed_emb = Ring.embed(normed).to_tuple()
        self.assertEqual((round(normed_emb[0], 3), round(normed_emb[1], 3)),
                         (round(best[0], 3), round(best[1], 3)),
                         msg=f"Standard imaginary normalization did not yield minimal embedding for {z}")

    
    @given(st.integers(min_value=2, max_value=50))
    def test_factorization_standard(self, n):
        x = self.Ring(n, 0)
        factors = self.Ring.factorize(x)
        prod = self.Ring.mul_id()
        for fac, exp in factors:
            for _ in range(int(exp)):
                prod = prod * fac
        self.assertEqual(prod.normalize(), x.normalize(),
                         msg=f"Standard: Factorization reconstruction failed for {x}")
    
    def test_identity_standard(self):
        self.assertEqual(self.Ring.zero(), self.Ring(0, 0),
                         msg="Standard: zero() incorrect in QuadIntRing")
        self.assertEqual(self.Ring.mul_id(), self.Ring(1, 0),
                         msg="Standard: mul_id() incorrect in QuadIntRing")

class TestQuadIntRingHalfIntegral(unittest.TestCase):
    """
    Tests for the factory QuadIntRing in the half-integral representation (e.g. d = 5).
    We test extended operations such as gcd, factorization, and normalization.
    """
    def setUp(self):
        self.d = 5
        self.Ring = QuadIntRing(self.d, force_ufd=True)
    
    @given(st.integers(min_value=-100, max_value=100))
    def test_gcd_half_integral(self, a):
        x = self.Ring(a, 0)
        y = self.Ring(a + 14, 0)
        g = self.Ring.gcd(x, y)
        expected = math.gcd(a, a + 14)
        self.assertEqual(abs(int(g.a)), expected,
                         msg=f"Half-integral: gcd mismatch for {a} and {a+14}")
        
    def test_fundamental_unit_half_integral(self):
        """
        For d = 5, the ring is half‑integral and the canonical fundamental unit is (1+√5)/2.
        Its real embedding (via the half‑integral transformation) should have absolute value
        approximately (1+√5)/2.
        """
        fu = self.Ring.get_fundamental_unit()
        emb = abs(float(self.Ring.embed(fu)))
        expected = (1 + math.sqrt(5)) / 2.0  # approximately 1.61803...
        self.assertAlmostEqual(emb, expected, delta=0.001,
            msg=f"Half-integral: Fundamental unit embedding mismatch for d=5: got {emb}, expected {expected}")
    
    def test_imaginary_units_minus3(self):
        d = -3
        Ring = QuadIntRing(d, force_ufd=True)
        units = Ring.get_imaginary_units()
        embeddings = {complex(*Ring.embed(u).to_tuple()) for u in units}
        # The canonical units in Z[√-3] are {1, -1, ω, -ω, ω^2, -ω^2},
        # where ω = (1+√-3)/2, which numerically equals:
        # ω = 0.5 + i*(√3/2)  and ω^2 = -0.5 + i*(√3/2).
        expected = {
            1+0j,
            -1+0j,
            0.5+1j*(math.sqrt(3)/2),
            -0.5-1j*(math.sqrt(3)/2),
            -0.5+1j*(math.sqrt(3)/2),
            0.5-1j*(math.sqrt(3)/2)
        }
        # Round to 3 decimal places to account for floating-point imprecision.
        rounded_embeddings = { (round(c.real, 3), round(c.imag, 3)) for c in embeddings }
        rounded_expected = { (round(c.real, 3), round(c.imag, 3)) for c in expected }
        self.assertEqual(rounded_embeddings, rounded_expected,
            msg=f"Imaginary: Units for d=-3 are incorrect: got {rounded_embeddings}")

        
    @given(st.integers(min_value=-50, max_value=50), st.integers(min_value=-50, max_value=50))
    def test_normalize_half_integral(self, a, b):
        assume(not (a == 0 and b == 0))
        z = self.Ring(a, b)
        zn = z.normalize()
        self.assertEqual(zn.normalize(), zn,
                         msg=f"Half-integral normalization not idempotent for element {z}")
        
    @given(st.integers(min_value=-50, max_value=50), st.integers(min_value=-50, max_value=50))
    def test_normalize_embedding_range_half_integral(self, a, b):
        assume(not (a == 0 and b == 0))
        z = self.Ring(a, b).normalize()
        if abs(int(QuadInt.norm(z))) == 1:
            return
        sigma_z = abs(float(self.Ring.embed(z)))
        epsilon = self.Ring.get_fundamental_unit()
        sigma_epsilon = abs(float(self.Ring.embed(epsilon)))
        self.assertGreaterEqual(sigma_z, 1,
                                msg=f"Half-integral: Embedding of normalized {z} is less than 1: {sigma_z}")
        self.assertLess(sigma_z, sigma_epsilon,
                        msg=f"Half-integral: Embedding of normalized {z} is not less than |σ(ε)| = {sigma_epsilon}")
        
    def test_normalize_idempotence_half_integral_imaginary(self):
        d = -3
        Ring = QuadIntRing(d, force_ufd=True)
        test_elements = [
            Ring(2, 3),    # Represents 2 + 3ω, where ω = (1+√(-3))/2
            Ring(-1, 4),
            Ring(5, -2)
        ]
        for z in test_elements:
            zn = z.normalize()
            self.assertEqual(zn.normalize(), zn,
                             msg=f"Half-integral imaginary normalization not idempotent for element {z}")

    def test_normalize_embedding_minimal_half_integral_imaginary(self):
        d = -3
        Ring = QuadIntRing(d, force_ufd=True)
        z = Ring(2, 3)  # Fixed test element.
        normed = z.normalize()
        U = Ring.get_imaginary_units()
        candidate_embeddings = [Ring.embed(u * z).to_tuple() for u in U]
        best = min(candidate_embeddings)
        normed_emb = Ring.embed(normed).to_tuple()
        self.assertEqual((round(normed_emb[0], 3), round(normed_emb[1], 3)),
                         (round(best[0], 3), round(best[1], 3)),
                         msg=f"Half-integral imaginary normalization did not yield minimal embedding for {z}")

    
    @given(st.integers(min_value=2, max_value=5000))
    def test_factorization_half_integral(self, n):
        x = self.Ring(n, 0)
        factors = self.Ring.factorize(x)
        prod = self.Ring.mul_id()
        for fac, exp in factors:
            for _ in range(int(exp)):
                prod = prod * fac
        self.assertEqual(prod.normalize(), x.normalize(),
                         msg=f"Half-integral: Factorization reconstruction failed for {x}")
    
    def test_identity_half_integral(self):
        self.assertEqual(self.Ring.zero(), self.Ring(0, 0),
                         msg="Half-integral: zero() incorrect in QuadIntRing")
        self.assertEqual(self.Ring.mul_id(), self.Ring(1, 0),
                         msg="Half-integral: mul_id() incorrect in QuadIntRing")

if __name__ == '__main__':
    unittest.main()