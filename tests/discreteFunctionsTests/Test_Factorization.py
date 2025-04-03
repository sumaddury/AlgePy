import unittest
from hypothesis import given, settings, strategies as st
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from DiscreteFunctions import Factorization, PrimalityTesting

#python -m unittest tests/discreteFunctionsTests/Test_Factorization.py

class TestFactorizationDivisors(unittest.TestCase):
    @given(st.integers(min_value=1, max_value=10000))
    def test_divisors_correctness(self, n):
        divs = Factorization.divisors(n)
        self.assertIn(1, divs, msg=f"1 not in divisors of {n}")
        self.assertIn(n, divs, msg=f"{n} not in divisors of {n}")
        for d in divs:
            self.assertEqual(n % d, 0, msg=f"{d} does not divide {n}")

class TestFactorizationPrimeFactors(unittest.TestCase):
    @given(st.integers(min_value=2, max_value=10000))
    def test_prime_factors_are_prime(self, n):
        pf = Factorization.prime_factors(n)
        for p in pf:
            self.assertTrue(PrimalityTesting.eratosthenes_primality(p),
                            msg=f"Prime factor {p} of {n} failed primality test")

class TestFactorizationReconstruction(unittest.TestCase):
    @given(st.integers(min_value=2, max_value=10000))
    def test_factorize_reconstructs(self, n):
        factors = Factorization.factorize(n)
        product = 1
        for prime, exp in factors:
            product *= prime ** exp
        self.assertEqual(product, n,
                         msg=f"Reconstruction failed for {n}: got {product}")

if __name__ == '__main__':
    unittest.main()
