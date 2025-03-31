import unittest
from hypothesis import given, settings, strategies as st
import sys, os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from SingletonStructures import Z, R, Z_n, Z_mod_
from DiscreteFunctions import PrimalityTesting, Factorization, ArithmeticFunctions

# --- Helper: Use DiscreteFunctions' coprime for checking unit property ---
def is_unit_in_Zn(a: int, n: int) -> bool:
    return PrimalityTesting.coprime(a, n)

def mod_repr(a: int, n: int) -> int:
    return a % n

#python -m unittest tests/Test_SingletonStructures.py

class TestZ(unittest.TestCase):
    @given(st.integers(min_value=-1000, max_value=1000), st.integers(min_value=-1000, max_value=1000))
    def test_addition_commutativity(self, a, b):
        A = Z(a)
        B = Z(b)
        self.assertEqual(A + B, B + A,
                         msg=f"Addition not commutative for {a} and {b}")

    @given(st.integers(min_value=-1000, max_value=1000))
    def test_negation_identity(self, a):
        A = Z(a)
        self.assertEqual(A + (-A), Z.zero(),
                         msg=f"Negation failed for {a}")

    @given(st.integers(min_value=-1000, max_value=1000), st.integers(min_value=0, max_value=10))
    def test_power(self, a, exp):
        if exp < 0:
            return
        A = Z(a)
        result = A.__pow__(Z(exp))
        self.assertEqual(result, Z(a ** exp),
                         msg=f"Power failed for a={a} and exp={exp}")

class TestR(unittest.TestCase):
    @given(st.floats(min_value=-1000, max_value=1000, allow_nan=False, allow_infinity=False),
           st.floats(min_value=-1000, max_value=1000, allow_nan=False, allow_infinity=False))
    def test_addition_commutativity(self, a, b):
        A = R(a)
        B = R(b)
        self.assertAlmostEqual(float(A + B), float(B + A), places=7,
                               msg=f"Addition not commutative for {a} and {b}")

    @given(st.floats(min_value=-1000, max_value=1000, allow_nan=False, allow_infinity=False))
    def test_division_inverse(self, a):
        if a == 0:
            return
        A = R(a)
        self.assertAlmostEqual(float(A / A), 1.0, places=7,
                               msg=f"Division inverse failed for {a}")

class TestZ_n(unittest.TestCase):
    @given(st.integers(min_value=1, max_value=100), st.integers(min_value=-1000, max_value=1000))
    def test_modular_reduction(self, n, a):
        A = Z_n(a, n)
        expected = mod_repr(a, n)
        self.assertEqual(int(A.z), expected,
                         msg=f"Modular reduction failed: {a} mod {n} != {expected}")

    @given(st.integers(min_value=2, max_value=100), st.integers(min_value=-1000, max_value=1000))
    def test_addition_mod(self, n, a):
        b = a + 17  # arbitrary offset
        A = Z_n(a, n)
        B = Z_n(b, n)
        result = A + B
        expected = mod_repr(a + b, n)
        self.assertEqual(int(result.z), expected,
                         msg=f"Modular addition failed for a={a}, b={b} mod {n}")

    @given(st.integers(min_value=2, max_value=100), st.integers(min_value=1, max_value=1000))
    def test_multiplicative_inverse(self, n, a):
        if not is_unit_in_Zn(a, n):
            return
        A = Z_n(a, n)
        invA = A.inverse()
        prod = A * invA
        self.assertEqual(int(prod.z), 1,
                         msg=f"Multiplicative inverse failed for a={a} mod {n}")

    @given(st.integers(min_value=2, max_value=100), st.integers(min_value=1, max_value=1000))
    def test_order_property(self, n, a):
        if not is_unit_in_Zn(a, n):
            return
        A = Z_n(a, n)
        ordA = A.order()
        result = A.__pow__(ordA)
        self.assertEqual(int(result.z), 1,
                         msg=f"Order computation failed for a={a} mod {n}, order {ordA}")

class TestZ_mod_Function(unittest.TestCase):
    def test_constructor_for_prime_modulus(self):
        prime_modulus = 7
        Zp = Z_mod_(prime_modulus)
        self.assertTrue(Zp.is_cyclic(),
                        msg=f"Expected Z_ for prime {prime_modulus} to be cyclic")
        for a in range(1, prime_modulus):
            elem = Zp(a)
            # Check that the unit condition (via PrimalityTesting.coprime) holds
            self.assertTrue(is_unit_in_Zn(a, prime_modulus),
                            msg=f"Expected {a} to be a unit mod {prime_modulus}")
            inv = elem.inverse()
            prod = elem * inv
            self.assertEqual(int(prod.z), 1,
                             msg=f"Inverse failed for a={a} mod {prime_modulus}")
        # Test legendre: For 3 mod 7, the Legendre symbol should be either -1, 0, or 1.
        result = Zp.legendre(Z(3))
        self.assertIn(int(result), [-1, 0, 1],
                      msg="Legendre symbol should be -1, 0, or 1")
    
    def test_constructor_for_composite_modulus(self):
        composite_modulus = 8
        Z8 = Z_mod_(composite_modulus)
        units = Z8.units()
        computed = sorted(int(u.z) for u in units)
        expected = sorted(a for a in range(1, composite_modulus) if is_unit_in_Zn(a, composite_modulus))
        self.assertEqual(computed, expected,
                         msg=f"Units of Z_8 incorrect: expected {expected}, got {computed}")

    @given(st.integers(min_value=2, max_value=100))
    def test_Z_mod_constructor_range(self, n):
        Zmod = Z_mod_(n)
        for elem in Zmod.elements():
            self.assertGreaterEqual(int(elem.z), 0,
                                    msg=f"Element {elem} is less than 0 in Z_mod_(n)")
            self.assertLess(int(elem.z), n,
                            msg=f"Element {elem} is not less than modulus {n}")

if __name__ == '__main__':
    unittest.main()
