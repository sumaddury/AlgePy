import unittest
from hypothesis import given, settings, strategies as st
import math
import random
import sys
import os

# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from algepy.SingletonStructures import Z, Z_n
from algepy.DiscreteFunctions import PrimalityTesting

#python -m unittest tests/singletonStructuresTests/Test_Z_n.py


def is_unit_in_Zn(a: int, n: int) -> bool:
    """
    Returns True if gcd(a,n) == 1.
    """
    return PrimalityTesting.coprime(a, n)

class TestZ_n(unittest.TestCase):

    @given(st.integers(min_value=1, max_value=200), st.integers(min_value=-1000, max_value=1000))
    def test_modular_reduction(self, n, a):
        """
        (a mod n) should match int(Z_n(a,n).z).
        """
        A = Z_n(a, n)
        self.assertEqual(int(A.z), a % n,
                         msg=f"Modular reduction failed: {a} mod {n}")

    @given(st.integers(min_value=2, max_value=100),
           st.integers(min_value=-100, max_value=100),
           st.integers(min_value=-100, max_value=100))
    def test_addition_mod(self, n, a, b):
        """
        (a+b) mod n = (int(Z_n(a,n).z) + int(Z_n(b,n).z)) mod n
        """
        A = Z_n(a, n)
        B = Z_n(b, n)
        R = A + B
        expected = ((a % n) + (b % n)) % n
        self.assertEqual(int(R.z), expected,
                         msg=f"Modular addition mismatch for {a}+{b} mod {n}")

    @given(st.integers(min_value=2, max_value=100),
           st.integers(min_value=-100, max_value=100),
           st.integers(min_value=-100, max_value=100))
    def test_multiplication_mod(self, n, a, b):
        A = Z_n(a, n)
        B = Z_n(b, n)
        R = A * B
        expected = ((a % n) * (b % n)) % n
        self.assertEqual(int(R.z), expected,
                         msg=f"Modular multiplication mismatch for {a}*{b} mod {n}")

    @given(st.integers(min_value=2, max_value=100),
           st.integers(min_value=1, max_value=2000))
    def test_inverse(self, n, a):
        """
        If gcd(a,n)=1, then (Z_n(a)*Z_n(a).inverse()).z == 1.
        """
        if not is_unit_in_Zn(a, n):
            return
        A = Z_n(a, n)
        invA = A.inverse()
        prod = A * invA
        self.assertEqual(int(prod.z), 1,
                         msg=f"Inverse failed for {a} mod {n}")

    @given(st.integers(min_value=2, max_value=50),
           st.integers(min_value=1, max_value=1000))
    def test_order(self, n, a):
        """
        order(a mod n) is the smallest k such that (a^k)=1 in Z_n
        """
        if not is_unit_in_Zn(a, n):
            return
        A = Z_n(a, n)
        k = A.order()
        # Check that a^k = 1 mod n
        self.assertEqual(int((A ** k).z), 1,
                         msg=f"Order computation failed for {a} mod {n}")

if __name__ == "__main__":
    unittest.main()
