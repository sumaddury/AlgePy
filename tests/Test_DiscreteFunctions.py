import unittest
from hypothesis import given, settings, strategies as st
import math
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from DiscreteFunctions import (
    PrimalityTesting,
    PrimeNumberTheorem,
    Factorization,
    ArithmeticFunctions
)

#python -m unittest tests/Test_DiscreteFunctions.py

# --- Helper: Known small primes (for deterministic tests) ---
KNOWN_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
# Known small non-Carmichael composites.
KNOWN_COMPOSITES = [4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25]

class TestPrimalityTesting(unittest.TestCase):
    @given(st.integers(min_value=1, max_value=1000), st.integers(min_value=1, max_value=1000))
    def test_coprime_property(self, a, b):
        self.assertEqual(PrimalityTesting.coprime(a, b), math.gcd(a, b) == 1,
                         msg=f"Failed for a={a}, b={b}")

    @given(st.integers(min_value=1, max_value=1000), st.integers(min_value=1, max_value=1000))
    def test_extended_euclidean_property(self, a, b):
        if a == 0 and b == 0:
            return
        x, y = PrimalityTesting.extended_euclidean(a, b)
        self.assertEqual(a * x + b * y, math.gcd(a, b),
                         msg=f"Failed for a={a}, b={b}: got x={x}, y={y}")

    def test_eratosthenes_on_known_numbers(self):
        for p in KNOWN_PRIMES:
            self.assertTrue(PrimalityTesting.eratosthenes_primality(p),
                            msg=f"Expected prime {p} to pass sieve test")
        for c in KNOWN_COMPOSITES:
            self.assertFalse(PrimalityTesting.eratosthenes_primality(c),
                             msg=f"Expected composite {c} to fail sieve test")

    @settings(max_examples=20)
    @given(st.sampled_from(KNOWN_PRIMES))
    def test_fermat_on_primes(self, p):
        self.assertTrue(PrimalityTesting.fermat_primality(p),
                        msg=f"Fermat test failed for prime {p}")

    @settings(max_examples=20)
    @given(st.sampled_from(KNOWN_COMPOSITES))
    def test_fermat_on_composites(self, c):
        # Repeat the test 20 times for the same composite; count failures.
        false_count = 0
        trials = 20
        for _ in range(trials):
            if not PrimalityTesting.fermat_primality(c):
                false_count += 1
        # For non-Carmichael composites, we expect most trials to yield False.
        self.assertGreater(false_count, 0.9 * trials,
                           msg=f"Fermat test passed too often for composite {c}: only {false_count}/{trials} failures")

    @settings(max_examples=20)
    @given(st.sampled_from(KNOWN_PRIMES))
    def test_miller_rabin_on_primes(self, p):
        self.assertTrue(PrimalityTesting.miller_rabin_primality(p),
                        msg=f"Miller-Rabin test failed for prime {p}")

    @settings(max_examples=20)
    @given(st.sampled_from(KNOWN_COMPOSITES))
    def test_miller_rabin_on_composites(self, c):
        self.assertFalse(PrimalityTesting.miller_rabin_primality(c),
                         msg=f"Miller-Rabin test incorrectly passed for composite {c}")
        
class TestPrimeNumberTheorem(unittest.TestCase):
    @given(st.integers(min_value=1000, max_value=20000))
    def test_pi_monotonicity(self, x):
        pi_x = PrimeNumberTheorem.π(x)
        pi_x_next = PrimeNumberTheorem.π(x+100)
        self.assertLessEqual(pi_x, pi_x_next*1.05,
                             msg=f"π({x})={pi_x} > π({x+100})={pi_x_next}")

    @given(st.integers(min_value=1000, max_value=20000))
    def test_approximation_accuracy(self, x):
        pi_x = PrimeNumberTheorem.π(x)
        approx = PrimeNumberTheorem.approx_pnt_prime_count(x)
        self.assertAlmostEqual(approx, pi_x, delta=0.25 * pi_x,
                               msg=f"For x={x}, expected π(x)≈{pi_x} but got {approx}")
        
class TestFactorization(unittest.TestCase):
    @given(st.integers(min_value=2, max_value=10000))
    def test_divisor_validity(self, n):
        divs = Factorization.divisors(n)
        self.assertIn(1, divs, msg=f"1 not in divisors of {n}")
        self.assertIn(n, divs, msg=f"{n} not in its own divisors")
        for d in divs:
            self.assertEqual(n % d, 0, msg=f"{d} is not a divisor of {n}")

    @given(st.integers(min_value=2, max_value=10000))
    def test_factorize_reconstruction(self, n):
        factors = Factorization.factorize(n)
        prod = 1
        for prime, exp in factors:
            prod *= prime ** exp
        self.assertEqual(prod, n,
                         msg=f"Failed reconstruction for {n}: got {prod}")
        
class TestArithmeticFunctions(unittest.TestCase):
    @given(st.integers(min_value=2, max_value=1000))
    def test_sigma_for_prime(self, n):
        if PrimalityTesting.eratosthenes_primality(n):
            sigma_val = ArithmeticFunctions.σ(n)
            self.assertEqual(sigma_val, n + 1,
                             msg=f"For prime {n}, expected σ(n) = {n+1} but got {sigma_val}")

    @given(st.integers(min_value=1, max_value=1000))
    def test_phi_basic(self, n):
        phi_val = ArithmeticFunctions.ϕ(n)
        if n == 1:
            self.assertEqual(phi_val, 1, msg="ϕ(1) should be 1")
        else:
            self.assertLess(phi_val, n,
                            msg=f"For n={n}, expected ϕ(n) < {n} but got {phi_val}")

    def test_I_and_N_functions(self):
        for n in range(1, 101):
            self.assertEqual(ArithmeticFunctions.I(n), 1 if n == 1 else 0,
                             msg=f"I({n}) expected {1 if n==1 else 0}")
            self.assertEqual(ArithmeticFunctions.N(n), n,
                             msg=f"N({n}) expected {n}")

    @given(st.integers(min_value=1, max_value=1000))
    def test_lambda_property(self, n):
        lambda_val = ArithmeticFunctions.λ(n)
        Omega_val = ArithmeticFunctions.Ω(n)
        expected = (-1) ** Omega_val
        self.assertEqual(lambda_val, expected,
                         msg=f"For n={n}, expected λ(n)= {expected}, got {lambda_val}")

    @given(st.integers(min_value=1, max_value=1000))
    def test_mu_property(self, n):
        mu_val = ArithmeticFunctions.μ(n)
        factors = Factorization.factorize(n)
        if all(exp == 1 for _, exp in factors):
            expected = (-1) ** len(factors)
        else:
            expected = 0
        self.assertEqual(mu_val, expected,
                         msg=f"For n={n}, expected μ(n)= {expected}, got {mu_val}")

    def test_sigma_k_known_values(self):
        test_cases = [
            (6, 0, 4),   # Divisors of 6: [1,2,3,6] so 1+1+1+1 = 4
            (6, 1, 12),  # 1+2+3+6 = 12
            (6, 2, 50),  # 1+4+9+36 = 50
            (7, 0, 2),   # 7 is prime: divisors [1,7] so 1+1 = 2
            (7, 1, 8),   # 1+7 = 8
            (7, 2, 50),  # 1+49 = 50
        ]
        for n, k, expected in test_cases:
            with self.subTest(n=n, k=k):
                result = ArithmeticFunctions.σ_k(n, k)
                self.assertEqual(result, expected,
                                 msg=f"σₖ({n},{k}) expected {expected}, got {result}")

    @given(st.integers(min_value=1, max_value=1000))
    def test_sigma_zero_property(self, n):
        sigma0 = ArithmeticFunctions.σ_k(n, 0)
        T_val = ArithmeticFunctions.Τ(n)
        self.assertEqual(sigma0, T_val,
                         msg=f"For n={n} with k=0, expected σ₀(n)=T(n)={T_val}, got {sigma0}")

    @given(st.integers(min_value=1, max_value=1000))
    def test_sigma_one_property(self, n):
        sigma1 = ArithmeticFunctions.σ_k(n, 1)
        sigma_val = ArithmeticFunctions.σ(n)
        self.assertEqual(sigma1, sigma_val,
                         msg=f"For n={n} with k=1, expected σ₁(n)=σ(n)={sigma_val}, got {sigma1}")

if __name__ == '__main__':
    unittest.main()