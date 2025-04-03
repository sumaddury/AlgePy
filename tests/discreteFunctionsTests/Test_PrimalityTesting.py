import unittest
from hypothesis import given, settings, strategies as st
import math
import random
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))
from DiscreteFunctions import PrimalityTesting

#python -m unittest tests/discreteFunctionsTests/Test_PrimalityTesting.py

CARMICHAEL_NUMS = {561, 1105, 1729, 2465, 2821, 6601, 8911}

def is_prime_det(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return n == 2
    limit = math.isqrt(n)
    for i in range(3, limit + 1, 2):
        if n % i == 0:
            return False
    return True

class TestCoprimeExtendedEuclidean(unittest.TestCase):
    @given(st.integers(min_value=1, max_value=10000), st.integers(min_value=1, max_value=10000))
    def test_coprime_vs_gcd(self, a, b):
        """
        Property: PrimalityTesting.coprime(a,b) <=> gcd(a,b) == 1
        """
        gcd_val = math.gcd(a, b)
        result = PrimalityTesting.coprime(a, b)
        self.assertEqual(result, (gcd_val == 1),
                         msg=f"coprime({a},{b}) mismatch with gcd={gcd_val}")

    @given(st.integers(min_value=0, max_value=10000), st.integers(min_value=0, max_value=10000))
    def test_extended_euclidean(self, a, b):
        """
        Property: a*x + b*y = gcd(a,b) for the returned x,y (except a=b=0).
        """
        if a == 0 and b == 0:
            return
        x, y = PrimalityTesting.extended_euclidean(a, b)
        gcd_val = math.gcd(a, b)
        lhs = a * x + b * y
        self.assertEqual(lhs, gcd_val,
                         msg=f"extended_euclidean({a},{b}) -> x={x}, y={y}, but a*x+b*y={lhs} != gcd={gcd_val}")

class TestEratosthenes(unittest.TestCase):

    @given(st.integers(min_value=2, max_value=10000))
    def test_eratosthenes_correctness(self, n):
        """
        Compare eratosthenes_primality(n) with our reference is_prime_det(n).
        """
        expected = is_prime_det(n)
        result = PrimalityTesting.eratosthenes_primality(n)
        self.assertEqual(result, expected,
                         msg=f"eratosthenes_primality({n}) mismatch. Expected {expected}.")

class TestFermatPrimality(unittest.TestCase):
    @settings(max_examples=10000)
    @given(st.lists(st.integers(min_value=2, max_value=100000), min_size=100, max_size=100))
    def test_fermat_aggregate_accuracy(self, numbers):
        """
        We collect a batch of numbers, separate primes and composites (excluding Carmichael).
        Then we measure how many times Fermat's test is correct.
        """
        primes = [n for n in numbers if is_prime_det(n)]
        comps = [n for n in numbers if not is_prime_det(n) and n not in CARMICHAEL_NUMS]

        trials = 5
        total = 0
        correct = 0

        for p in primes:
            for _ in range(trials):
                total += 1
                if PrimalityTesting.fermat_primality(p):
                    correct += 1

        for c in comps:
            for _ in range(trials):
                total += 1
                if not PrimalityTesting.fermat_primality(c):
                    correct += 1

        if total == 0:
            return
        accuracy = correct / total
        self.assertGreaterEqual(accuracy, 0.98,
            msg=f"Fermat overall accuracy too low: {accuracy:.2f} over {total} checks.")

class TestMillerRabinPrimality(unittest.TestCase):

    @settings(max_examples=10000)
    @given(st.lists(st.integers(min_value=2, max_value=100000), min_size=100, max_size=100))
    def test_miller_rabin_aggregate_accuracy(self, numbers):
        """
        We collect a batch of numbers, separate primes and composites (excluding Carmichael).
        Then measure how many times Miller-Rabin is correct.
        """
        primes = [n for n in numbers if is_prime_det(n)]
        comps = [n for n in numbers if not is_prime_det(n) and n not in CARMICHAEL_NUMS]

        trials = 5
        total = 0
        correct = 0

        for p in primes:
            for _ in range(trials):
                total += 1
                if PrimalityTesting.miller_rabin_primality(p):
                    correct += 1

        for c in comps:
            for _ in range(trials):
                total += 1
                if not PrimalityTesting.miller_rabin_primality(c):
                    correct += 1

        if total == 0:
            return
        accuracy = correct / total
        self.assertGreaterEqual(accuracy, 0.98,
            msg=f"Miller-Rabin overall accuracy too low: {accuracy:.2f} over {total} checks.")

if __name__ == '__main__':
    unittest.main()
