import unittest
from hypothesis import given, assume, settings, strategies as st, HealthCheck
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from .DiscreteFunctions import ArithmeticFunctions, Factorization, PrimalityTesting

#python -m unittest tests/discreteFunctionsTests/Test_ArithmeticFunctions.py

class TestSigmaFunctions(unittest.TestCase):
    @given(st.integers(min_value=1, max_value=10000))
    def test_sigma_definition(self, n):
        divs = Factorization.divisors(n)
        self.assertEqual(ArithmeticFunctions.σ(n), sum(divs),
                         msg=f"σ({n}) mismatch.")

    @given(st.integers(min_value=1, max_value=10000), st.integers(min_value=-3, max_value=3))
    def test_sigma_k_consistency(self, n, k):
        if k == 0:
            self.assertEqual(ArithmeticFunctions.σ_k(n, 0), len(Factorization.divisors(n)),
                             msg=f"σ₀({n}) mismatch.")
        elif k == 1:
            self.assertEqual(ArithmeticFunctions.σ_k(n, 1), ArithmeticFunctions.σ(n),
                             msg=f"σ₁({n}) mismatch.")

class TestTotientFunctions(unittest.TestCase):
    @settings(suppress_health_check=[HealthCheck.filter_too_much])
    @given(st.integers(min_value=2, max_value=10000))
    def test_totient_prime(self, n):
        assume(PrimalityTesting.eratosthenes_primality(n))
        self.assertEqual(ArithmeticFunctions.ϕ(n), n - 1,
                         msg=f"ϕ({n}) should be {n-1} for a prime, got {ArithmeticFunctions.ϕ(n)}")

    @given(st.integers(min_value=1, max_value=10000), st.integers(min_value=1, max_value=1000))
    def test_totient_multiplicativity(self, m, n):
        assume(ArithmeticFunctions.ϕ(m) > 0 and ArithmeticFunctions.ϕ(n) > 0)
        assume(math.gcd(m, n) == 1)
        self.assertEqual(ArithmeticFunctions.ϕ(m * n), ArithmeticFunctions.ϕ(m) * ArithmeticFunctions.ϕ(n),
                         msg=f"ϕ multiplicativity failed for {m} and {n}")

class TestMöbiusAndLiouville(unittest.TestCase):
    @given(st.integers(min_value=1, max_value=10000))
    def test_mu_zero_when_not_square_free(self, n):
        factors = Factorization.factorize(n)
        if any(exp > 1 for _, exp in factors):
            self.assertEqual(ArithmeticFunctions.μ(n), 0,
                             msg=f"μ({n}) expected 0 for non-square-free, got {ArithmeticFunctions.μ(n)}")

    @given(st.integers(min_value=1, max_value=10000))
    def test_mu_values_for_square_free(self, n):
        factors = Factorization.factorize(n)
        if all(exp == 1 for _, exp in factors):
            expected = (-1) ** len(factors)
            self.assertEqual(ArithmeticFunctions.μ(n), expected,
                             msg=f"μ({n}) mismatch for square-free number.")

    @given(st.integers(min_value=1, max_value=10000))
    def test_liouville_function(self, n):
        omega = ArithmeticFunctions.Ω(n)
        expected = (-1) ** omega
        self.assertEqual(ArithmeticFunctions.λ(n), expected,
                         msg=f"λ({n}) expected {expected} but got {ArithmeticFunctions.λ(n)}")

class TestIdentityFunctions(unittest.TestCase):
    @given(st.integers(min_value=1, max_value=10000))
    def test_I(self, n):
        self.assertEqual(ArithmeticFunctions.I(n), 1 if n == 1 else 0,
                         msg=f"I({n}) incorrect.")
    
    @given(st.integers(min_value=1, max_value=10000))
    def test_N(self, n):
        self.assertEqual(ArithmeticFunctions.N(n), n,
                         msg=f"N({n}) incorrect.")
        
class TestArithmeticFunctions_PerfectSquareFree(unittest.TestCase):
    def test_is_perfect_known_values(self):
        known_perfect = [6, 28, 496]
        for n in known_perfect:
            self.assertTrue(ArithmeticFunctions.is_perfect(n),
                            msg=f"{n} should be perfect but is_perfect({n}) returned False.")
        
        non_perfect = [8, 12, 18, 20, 24, 30, 100, 200]
        for n in non_perfect:
            self.assertFalse(ArithmeticFunctions.is_perfect(n),
                             msg=f"{n} is not perfect but is_perfect({n}) returned True.")

    def test_is_square_free_known_values(self):
        square_free = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29]
        for n in square_free:
            self.assertTrue(ArithmeticFunctions.is_square_free(n),
                            msg=f"{n} is square-free but is_square_free({n}) returned False.")
        
        non_square_free = [4, 8, 9, 12, 16, 18, 20, 24, 25, 27, 28, 32, 36, 49]
        for n in non_square_free:
            self.assertFalse(ArithmeticFunctions.is_square_free(n),
                             msg=f"{n} is not square-free but is_square_free({n}) returned True.")

if __name__ == '__main__':
    unittest.main()
