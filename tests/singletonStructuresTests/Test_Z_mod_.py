import unittest
from hypothesis import given, strategies as st
import math
import random
import sys
import os

# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'algepy')))

from algepy.SingletonStructures import Z_mod_
from algepy.DiscreteFunctions import PrimalityTesting, ArithmeticFunctions

#python -m unittest tests/singletonStructuresTests/Test_Z_mod_.py

class TestZ_modConstructor(unittest.TestCase):

    def test_prime_modulus(self):
        prime_mod = 7
        Z7 = Z_mod_(prime_mod)
        self.assertTrue(Z7.is_cyclic(),
                        msg="Expected Z_7 to be cyclic for prime modulus=7.")
        self.assertEqual(Z7.unit_order().z, prime_mod - 1,
                         msg="unit_order mismatch for prime modulus 7.")

        if prime_mod > 2:
            symbol = Z7.legendre(Z7(3))
            self.assertIn(symbol.z, {-1, 0, 1},
                          msg="Legendre symbol should be -1, 0, or 1")

    def test_composite_modulus(self):
        composite_mod = 12
        Z12 = Z_mod_(composite_mod)

        self.assertIsInstance(Z12.is_cyclic(), bool,
                              msg="is_cyclic() for composite should return a bool")

        self.assertEqual(Z12.unit_order().z, 4,
                         msg=f"unit_order mismatch for composite modulus {composite_mod}.")

        all_elems = Z12.elements()
        self.assertEqual(len(all_elems), composite_mod,
                         msg=f"elements() mismatch for modulus {composite_mod}.")
        all_units = Z12.units()
        self.assertEqual(len(all_units), 4,
                         msg=f"units() mismatch for modulus {composite_mod} (ϕ(12)=4).")

    @given(st.integers(min_value=2, max_value=50))
    def test_range_of_moduli(self, n):
        """
        For a range of moduli, ensure Z_mod_(n).elements() has length n,
        and Z_mod_(n).units() matches ϕ(n).
        """
        if n == 0:
            return
        ZN = Z_mod_(n)
        elems = ZN.elements()
        self.assertEqual(len(elems), n,
                         msg=f"Expected {n} elements in Z_{n}, got {len(elems)}")

        # Compare len(units) to ArithmeticFunctions.ϕ(n) if n>1
        if n > 1:
            phi_n = ArithmeticFunctions.ϕ(n)
            self.assertEqual(len(ZN.units()), phi_n,
                             msg=f"units() mismatch for modulus {n}")

if __name__ == "__main__":
    unittest.main()
