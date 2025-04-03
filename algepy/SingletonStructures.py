import typing
from typing import Callable, Self
import math
import sys
from .DiscreteFunctions import PrimalityTesting, Factorization, ArithmeticFunctions

class Z:
    r"""
    Represents an element of the ring of integers, ℤ.

    This class wraps a Python integer and implements standard arithmetic operations 
    (addition, multiplication, subtraction, etc.) consistent with the algebraic structure of ℤ.

    **Preconditions:**
    
    - The input`z` must be a Python integer.
    
    **Mathematical Intuition:**
    
    The set ℤ is closed under addition, subtraction, and multiplication. In particular, for any 
    integers :math:`a, b \in \mathbb{Z}`:
    
    .. math::
        a + b \in \mathbb{Z}, \quad a - b \in \mathbb{Z}, \quad a \cdot b \in \mathbb{Z}.
    
    This class provides operator overloading so that arithmetic in ℤ can be performed using natural syntax.
    
    :param z: An integer representing the element in ℤ.
    :raises Exception: If`z` is not an integer.
    """
    def __init__(self, z: int) -> Self:
        if isinstance(z, int):
            self.z = z
        else:
            raise Exception("z ¬∈ Z")

    def __add__(self, b: Self) -> Self:
        return Z(self.z + b.z)

    def __mul__(self, b: Self) -> Self:
        return Z(self.z * b.z)

    def __neg__(self) -> Self:
        return Z(-1 * self.z)

    def __sub__(self, b: Self) -> Self:
        return Z(self.z - b.z)

    def __pow__(self, b: Self) -> Self:
        return Z(self.z ** b.z)

    def __floordiv__(self, b: Self) -> Self:
        return Z(self.z // b.z)

    def __mod__(self, b: Self) -> Self:
        return Z(self.z % b.z)

    def __eq__(self, b: Self) -> bool:
        return self.z == b.z

    def __repr__(self) -> str:
        return str(self.z)

    def __int__(self) -> int:
        return self.z

    @staticmethod
    def zero() -> Self:
        r"""
        Returns the additive identity (0) in ℤ.
        
        **Mathematical Intuition:**
        
        The identity element :math:`0 \in \mathbb{Z}` satisfies:
        
        .. math::
            a + 0 = a \quad \forall a \in \mathbb{Z}.
        
        :return: The integer 0 wrapped as a Z object.
        """
        return Z(0)

    @staticmethod
    def mul_id() -> Self:
        r"""
        Returns the multiplicative identity (1) in ℤ.
        
        **Mathematical Intuition:**
        
        The identity element :math:`1 \in \mathbb{Z}` satisfies:
        
        .. math::
            a \cdot 1 = a \quad \forall a \in \mathbb{Z}.
        
        :return: The integer 1 wrapped as a Z object.
        """
        return Z(1)

    @staticmethod
    def abs(z: Self) -> Self:
        r"""
        Returns the absolute value of the integer element.
        
        **Mathematical Intuition:**
        
        The absolute value function is defined by:
        
        .. math::
            |a| = \begin{cases} a, & a \ge 0 \\ -a, & a < 0 \end{cases}.
        
        :param z: A Z object.
        :return: A new Z object with the absolute value of z.
        """
        return Z(abs(int(z)))

class R:
    r"""
    Represents an element of the field of real numbers, ℝ.

    This class wraps a real number (float or int) and defines standard arithmetic operations 
    as they occur in ℝ.

    **Preconditions:**
    
    - The input`r` must be a real number (either an int or a float).
    
    **Mathematical Intuition:**
    
    ℝ is a field; every nonzero element has a multiplicative inverse. For any :math:`r, s \in \mathbb{R}`:
    
    .. math::
        r + s \in \mathbb{R}, \quad r \cdot s \in \mathbb{R}, \quad r - s \in \mathbb{R}, \quad \text{and if } r \neq 0,\; \frac{1}{r} \in \mathbb{R}.
    
    This class models ℝ using floating-point arithmetic.
    
    :param r: A real number to be encapsulated as an element of ℝ.
    :raises Exception: If`r` is not a number.
    """
    def __init__(self, r: float) -> Self:
        if isinstance(r, float) or isinstance(r, int):
            self.r = float(r)
        else:
            raise Exception("r ¬∈ R")

    def __add__(self, b: Self) -> Self:
        return R(self.r + b.r)

    def __mul__(self, b: Self) -> Self:
        return R(self.r * b.r)

    def __neg__(self) -> Self:
        return R(-1 * self.r)

    def __sub__(self, b: Self) -> Self:
        return R(self.r - b.r)

    def __truediv__(self, b: Self) -> Self:
        return R(self.r / b.r)

    def __pow__(self, b: Self) -> Self:
        return R(self.r ** b.r)

    def __eq__(self, b: Self) -> bool:
        return self.r == b.r

    def __repr__(self) -> str:
        return str(self.r)

    def __float__(self) -> float:
        return self.r

    @staticmethod
    def zero() -> Self:
        r"""
        Returns the additive identity (0.0) in ℝ.
        
        **Mathematical Intuition:**
        
        0 is the unique element satisfying:
        
        .. math::
            r + 0 = r \quad \forall r \in \mathbb{R}.
        
        :return: An R object representing 0.
        """
        return R(0)

    @staticmethod
    def mul_id() -> Self:
        r"""
        Returns the multiplicative identity (1.0) in ℝ.
        
        **Mathematical Intuition:**
        
        1 is the unique element satisfying:
        
        .. math::
            r \cdot 1 = r \quad \forall r \in \mathbb{R}.
        
        :return: An R object representing 1.
        """
        return R(1)

    @staticmethod
    def abs(a: Self) -> Self:
        r"""
        Returns the absolute value of a real number element.
        
        **Mathematical Intuition:**
        
        Defined as:
        
        .. math::
            |r| = \begin{cases} r, & r \ge 0 \\ -r, & r < 0 \end{cases}.
        
        :param a: An R object.
        :return: An R object containing the absolute value of a.
        """
        return R(abs(float(a)))

class Z_n:
    r"""
    Represents an element of the ring of integers modulo n, ℤₙ.

    Each element is stored as a residue class modulo n. Arithmetic operations are performed 
    modulo n. This class supports addition, multiplication, subtraction, exponentiation, and, 
    for units, computing inverses and orders.

    **Preconditions:**
    
    - The modulus`n` must be a positive integer.
    - The element`z` must be a Python integer; it is automatically reduced modulo n.
    
    **Mathematical Intuition:**
    
    In ℤₙ, for an integer :math:`a` and modulus :math:`n`:
    
    .. math::
        a \mod n \in \{0, 1, \dots, n-1\},
    
    and when n is prime, ℤₙ forms a field. This class models the residue classes with the usual 
    modular arithmetic.
    
    :param z: An integer representing the element in ℤ.
    :param n: The modulus for the residue class.
    :raises Exception: If z or n is not an integer, or if n is invalid.
    """
    def __init__(self, z, n):
        try:
            self.z = z % n if isinstance(z, Z) else Z(z % int(n))
            self.n = n if isinstance(n, Z) else Z(n)
        except:
            raise Exception("z ¬∈ Z ∨ n ¬∈ Z")

    def __add__(self, b: Self) -> Self:
        if not self.n == b.n:
            raise Exception("unsupported operation")
        return self.__class__(self.z + b.z, self.n)

    def __mul__(self, b: Self) -> Self:
        if not self.n == b.n:
            raise Exception("unsupported operation")
        return self.__class__(self.z * b.z, self.n)

    def __neg__(self) -> Self:
        return self.__class__(-self.z, self.n)

    def __sub__(self, b: Self) -> Self:
        return self + -b

    def __pow__(self, b: Z) -> Self:
        return self.__class__(pow(int(self.z), int(b), int(self.n)),self.n)

    def __eq__(self, b: Self) -> bool:
        if not self.n == b.n:
            raise Exception("unsupported operation")
        return self.z == b.z

    def __repr__(self) -> str:
        return str(self.z) + " mod " + str(self.n)
    
    def order(self) -> Z:
        r"""
        Compute the order of the element in the multiplicative group of units U(n).

        **Mathematical Intuition:**
        
        For an element :math:`a \in \mathbb{Z}_n` that is invertible, its order is the smallest 
        positive integer :math:`k` such that:
        
        .. math::
            a^k \equiv 1 \pmod{n}.
        
        **Preconditions:**
        
        - The element must be a unit (i.e., :math:`\gcd(a, n) = 1`).
        
        :return: A Z object representing the order.
        :raises Exception: If the element is not a unit.
        """
        if not PrimalityTesting.coprime(int(self.z), int(self.n)):
            raise Exception(str(self)+" ¬∈ U_"+str(self.n))
        k = 1
        x = self
        while not x == Z_n(1,self.n):
            x *= self
            k += 1
        return Z(k)

    def inverse(self) -> Self:
        r"""
        Compute the multiplicative inverse of the element in ℤₙ.

        **Mathematical Intuition:**
        
        For a unit :math:`a \in \mathbb{Z}_n` (with :math:`\gcd(a, n) = 1`), there exists an inverse :math:`a^{-1}` 
        such that:
        
        .. math::
            a \cdot a^{-1} \equiv 1 \pmod{n}.
        
        **Preconditions:**
        
        - The element must be a unit.
        
        :return: A Z_n object representing the inverse.
        :raises Exception: If the element is not invertible.
        """
        if not PrimalityTesting.coprime(int(self.z), int(self.n)):
            raise Exception(str(self)+" ¬∈ U_"+str(self.n))
        return Z_n(PrimalityTesting.extended_euclidean(int(self.z), int(self.n))[0], self.n)

def Z_mod_(n: int) -> Callable[[int], Z_n]:
    r"""
    Returns a specialized constructor for the ring ℤₙ (integers modulo n).

    **Mathematical Intuition:**
    
    ℤₙ consists of residue classes:
    
    .. math::
        \mathbb{Z}_n = \{0, 1, \dots, n-1\},
    
    with operations defined modulo n. For prime n, ℤₙ is a field; for composite n, the unit group 
    (the set of invertible elements) is given by:
    
    .. math::
        U(n) = \{ a \in \mathbb{Z}_n : \gcd(a, n) = 1 \}.
    
    This function returns a subclass of Z_n that encapsulates additional properties (e.g., cyclicity, primitive roots)
    depending on whether n is prime or composite.
    
    **Preconditions:**
    
    - n must be a positive integer (n ∈ ℕ, n > 0).
    - n = 0 is not allowed.
    
    :param n: A positive integer representing the modulus.
    :return: A callable (class) for constructing elements of ℤₙ with the given modulus.
    :raises Exception: If n is not a positive integer or if n equals 0.
    """
    if not isinstance(n, int) or n < 1:
        raise Exception("n ¬∈ N")

    if n == 0:
        raise Exception("n = 0")
    
    class Z_(Z_n):
        """
        NOTE: This is the subclass constructed by function Z_mod_().
        """
        modulus = Z(n)

        def __init__(self, z):
            super().__init__(z, self.modulus)

        @classmethod
        def elements(cls) -> list[Self]:
            return [cls(k) for k in range(int(cls.modulus))]

        @classmethod
        def units(cls) -> list[Self]:
            return [cls(k) for k in range(1, int(cls.modulus)) if PrimalityTesting.coprime(k, int(cls.modulus))]

        @classmethod
        def zero(cls) -> Self:
            return cls(0)

        @classmethod
        def mul_id(cls) -> Self:
            return cls(1)

        @classmethod
        def total_order(cls) -> Z:
            return cls.modulus

        if not PrimalityTesting.eratosthenes_primality(n):
            @classmethod
            def unit_order(cls) -> Z:
                return Z(ArithmeticFunctions.ϕ(int(cls.modulus)))

            @classmethod
            def is_cyclic(cls) -> bool:
                r"""
                Determines whether ℤₙ is cyclic, i.e., whether its unit group is cyclic.

                **Mathematical Intuition:**
                
                A group G is cyclic if there exists a generator :math:`g` such that:
                
                .. math::
                    G = \{ g^k : k \in \mathbb{Z} \}.
                
                For ℤₙ, the cyclicity of the unit group depends on the prime factorization of n.
                
                **Preconditions:**
                
                - The modulus n must be greater than 1.
                
                :return: True if ℤₙ is cyclic, False otherwise.
                """
                if int(cls.modulus) in {1, 2, 4}:
                    return True
                factors = Factorization.prime_factors(int(cls.modulus))
                if len(factors) == 1 and factors[0] != 2 and PrimalityTesting.eratosthenes_primality(factors[0]):
                    return True
                if len(factors) == 2 and factors[0] == 2 and PrimalityTesting.eratosthenes_primality(factors[1]):
                    return True
                return False

            @classmethod
            def primitive_roots(cls) -> Self:
                r"""
                Finds the primitive roots of ℤₙ, if they exist.
                
                **Mathematical Intuition:**
                
                A primitive root modulo n is an element :math:`g \in U(n)` such that:
                
                .. math::
                    U(n) = \{ g^k : 0 \leq k < \varphi(n) \}.
                
                This function returns the list of all such generators when ℤₙ is cyclic.
                
                **Preconditions:**
                
                - ℤₙ must be cyclic.
                
                :return: A list of primitive roots in ℤₙ, or None if ℤₙ is not cyclic.
                :raises Exception: If ℤₙ is not cyclic.
                """
                if not cls.is_cyclic():
                    return None

                phi = cls.unit_order()
                factors = Factorization.prime_factors(int(phi))
                return [g for g in cls.units() if all(pow(int(g.z), int(phi) // q, int(cls.modulus)) != 1 for q in factors)]

        else:
            @classmethod
            def unit_order(cls) -> Z:
                return cls.modulus - Z(1)

            @classmethod
            def is_cyclic(cls) -> bool:
                return True

            @classmethod
            def primitive_roots(cls) -> Self:
                phi = cls.unit_order()
                factors = Factorization.prime_factors(int(phi))
                return [g for g in cls.units() if all(pow(int(g.z), int(phi) // q, int(cls.modulus)) != 1 for q in factors)]

            @classmethod
            def legendre(cls, z: Self) -> Z:
                r"""
                Computes the Legendre symbol :math:`\left(\frac{a}{p}\right)` for a in ℤ, where p is an odd prime.
                
                **Mathematical Intuition:**
                
                For an odd prime :math:`p`, the Legendre symbol is defined by:
                
                .. math::
                    \left(\frac{a}{p}\right) = \begin{cases}
                    0, & \text{if } p \mid a, \\
                    1, & \text{if } a \text{ is a quadratic residue modulo } p, \\
                    -1, & \text{if } a \text{ is a non-residue modulo } p.
                    \end{cases}
                
                **Preconditions:**
                
                - p (here, the modulus) must be an odd prime.
                
                :param z: An element of ℤₙ.
                :return: The Legendre symbol as a Z object.
                :raises Exception: If the modulus is 2.
                """
                if int(cls.modulus) == 2:
                    raise Exception("legendre symbol not defined for p=2")
                symbol = pow(int(z.z), (int(cls.modulus) - 1) // 2, int(cls.modulus))
                return Z(-1) if symbol == int(cls.modulus) - 1 else Z(symbol)

            def __truediv__(self, b: Self) -> Self:
                if not b in self.__class__.units():
                    raise Exception("unsupported operation")
                return self * b.inverse()

    Z_.__module__ = __name__
    if 'sphinx' in sys.modules:
        globals()[f'Z_'] = Z_

    return Z_

Z__ = Z_mod_(2)