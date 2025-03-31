import typing
from typing import Callable, Self
import math
from DiscreteFunctions import PrimalityTesting, Factorization, ArithmeticFunctions

class Z:
    """
    Represents an element of the ring of integers, ℤ.

    This class wraps a Python integer and implements standard arithmetic operations (addition, 
    multiplication, subtraction, etc.) consistent with the algebraic structure of ℤ.

    **Preconditions:**
    
    - The input `z` must be a Python integer.
    
    **Mathematical Intuition:**
    
    - The set ℤ is closed under addition, multiplication, and subtraction.
    - This class allows you to work with integers as algebraic objects with operator overloading 
      mimicking the standard operations in ℤ.
    
    :param z: An integer representing the element in ℤ.
    :raises Exception: If `z` is not an integer.
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
        """
        Returns the additive identity (0) in ℤ.
        
        :return: The integer 0 wrapped as a Z object.
        """
        return Z(0)

    @staticmethod
    def mul_id() -> Self:
        """
        Returns the multiplicative identity (1) in ℤ.
        
        :return: The integer 1 wrapped as a Z object.
        """
        return Z(1)

    @staticmethod
    def abs(z: Self) -> Self:
        """
        Returns the absolute value of the integer element.
        
        :param z: A Z object.
        :return: A new Z object with the absolute value of z.
        """
        return Z(abs(int(z)))

    @staticmethod
    def is_perfect(z: Self) -> bool:
        """
        Determines whether the integer is perfect.

        **Mathematical Intuition:**
        
        - A perfect number equals the sum of its proper divisors. Here we test whether the sum 
          of all divisors equals 2 * z.
        
        **Preconditions:**
        
        - The input must be a positive integer.
        
        :param z: A Z object.
        :return: True if z is a perfect number, otherwise False.
        """
        return ArithmeticFunctions.σ(int(z)) == 2 * int(z)

    @staticmethod
    def is_square_free(z: Self) -> bool:
        """
        Checks whether the integer is square-free (i.e., no prime factor appears with exponent > 1).

        **Preconditions:**
        
        - The input must be a positive integer.
        
        **Mathematical Intuition:**
        
        - A square-free number is not divisible by any perfect square greater than 1.
        
        :param z: A Z object.
        :return: True if z is square-free; False otherwise.
        """
        if int(z) < 1:
            return False
        if int(z) == 1:
            return True

        for p in range(2, math.isqrt(int(z))+1):
            if PrimalityTesting.eratosthenes_primality(p) and int(z) % (p*p) == 0:
                return False
        return True  

class R:
    """
    Represents an element of the field of real numbers, ℝ.

    This class wraps a real number (float or int) and defines standard arithmetic operations
    as they occur in ℝ.

    **Preconditions:**
    
    - The input `r` must be a real number (either an int or float).
    
    **Mathematical Intuition:**
    
    - ℝ is a field: every nonzero element has a multiplicative inverse, and all the standard 
      arithmetic operations are defined.
    
    :param r: A real number to be encapsulated as an element of ℝ.
    :raises Exception: If `r` is not a number.
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
        """
        Returns the additive identity (0.0) in ℝ.
        
        :return: An R object representing 0.
        """
        return R(0)

    @staticmethod
    def mul_id() -> Self:
        """
        Returns the multiplicative identity (1.0) in ℝ.
        
        :return: An R object representing 1.
        """
        return R(1)

    @staticmethod
    def abs(a: Self) -> Self:
        """
        Returns the absolute value of a real number element.
        
        :param a: An R object.
        :return: An R object containing the absolute value of a.
        """
        return R(abs(float(a)))

class Z_n:
    """
    Represents an element of the ring of integers modulo n, ℤₙ.

    Each element is stored as a residue class modulo n. Arithmetic operations are defined modulo n.
    This class supports addition, multiplication, subtraction, exponentiation, and, for units, 
    computing inverses and orders.

    **Preconditions:**
    
    - The modulus `n` must be a positive integer.
    - The element `z` must be an integer; it is automatically reduced modulo n.
    
    **Mathematical Intuition:**
    
    - ℤₙ forms a ring; when n is prime, ℤₙ forms a field. In ℤₙ, operations are performed modulo n.
    - The class includes methods to compute the order of an element (when it is invertible) 
      and to determine the multiplicative inverse.
    
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
        return Z_n(self.z + b.z, self.n)

    def __mul__(self, b: Self) -> Self:
        if not self.n == b.n:
            raise Exception("unsupported operation")
        return Z_n(self.z * b.z, self.n)

    def __neg__(self) -> Self:
        return Z_n(-self.z, self.n)

    def __sub__(self, b: Self) -> Self:
        return self + -b

    def __pow__(self, b: Z) -> Self:
        return Z_n(pow(int(self.z), int(b), int(self.n)),self.n)

    def __eq__(self, b: Self) -> bool:
        if not self.n == b.n:
            raise Exception("unsupported operation")
        return self.z == b.z

    def __repr__(self) -> str:
        return str(self.z) + " mod " + str(self.n)

    def order(self) -> Z:
        """
        Compute the order of the element in the multiplicative group of units U(n).

        **Preconditions:**
        
        - The element must be a unit in ℤₙ (i.e., coprime with n).
        
        **Mathematical Intuition:**
        
        - The order of an element is the smallest positive integer k such that the k-th power 
          of the element equals the multiplicative identity.
        
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
        """
        Compute the multiplicative inverse of the element in ℤₙ.

        **Preconditions:**
        
        - The element must be a unit in ℤₙ (i.e., coprime with n).
        
        :return: A Z_n object representing the multiplicative inverse.
        :raises Exception: If the element is not invertible.
        """
        if not PrimalityTesting.coprime(int(self.z), int(self.n)):
            raise Exception(str(self)+" ¬∈ U_"+str(self.n))
        return Z_n(PrimalityTesting.extended_euclidean(int(self.z), int(self.n))[0], self.n)

def Z_mod_(n: int) -> Callable[[int], Z_n]:
    """
    Returns a specialized constructor for the ring ℤₙ (integers modulo n).

    Depending on whether n is prime or composite, this function returns a subclass of Z_n with 
    additional class methods to list elements, units, compute orders, and (when applicable) determine 
    properties like cyclicity and primitive roots.

    NOTE: Read the source code for class methods the constructor provides. Sphinx does not document those.

    **Preconditions:**
    
    - n must be a positive integer (n ∈ ℕ and n > 0).
    - n = 0 is not allowed.
    
    **Mathematical Intuition:**
    
    - For composite moduli, ℤₙ may have a more complicated structure with non-trivial units.
    - For prime moduli, ℤₙ is a field where every nonzero element is invertible.
    - This function encapsulates the differences by returning an appropriate subclass of Z_n.
    
    :param n: A positive integer representing the modulus.
    :return: A callable (class) that constructs elements in ℤₙ with the given modulus.
    :raises Exception: If n is not a positive integer or if n equals 0.
    """
    if not isinstance(n, int) or n < 1:
        raise Exception("n ¬∈ N")

    if n == 0:
        raise Exception("n = 0")

    if not PrimalityTesting.eratosthenes_primality(n):
        class Z_(Z_n):
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

            @classmethod
            def unit_order(cls) -> Z:
                return Z(ArithmeticFunctions.ϕ(int(cls.modulus)))

            @classmethod
            def is_cyclic(cls) -> bool:
                """
                Determines whether ℤₙ is cyclic, i.e., whether it has a generator for its multiplicative group of units.

                **Mathematical Intuition:**
                
                - A group is cyclic if there exists an element (called a generator) such that every element of the group can be written as a power of this generator.
                - For ℤₙ, this property depends on the structure of the modulus n.
                
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
                """
                Finds the primitive roots of ℤₙ, if they exist.

                **Mathematical Intuition:**
                
                - A primitive root modulo n is an element of ℤₙ that generates all other units in ℤₙ under multiplication.
                - If ℤₙ is cyclic, the primitive roots are the elements that can produce all other units by exponentiation.
                
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
        class Z_(Z_n):
            modulus = Z(n)
    
            def __init__(self, z):
                super().__init__(z, self.modulus)
    
            @classmethod
            def elements(cls) -> list[Self]:
                return [cls(k) for k in range(int(cls.modulus))]

            @classmethod
            def units(cls) -> list[Self]:
                return cls.elements()[1:]

            @classmethod
            def zero(cls) -> Self:
                return cls(0)
    
            @classmethod
            def mul_id(cls) -> Self:
                return cls(1)

            @classmethod
            def total_order(cls) -> Z:
                return cls.modulus

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
                """
                Computes the Legendre symbol (a / p), which indicates whether a is a quadratic residue modulo p.

                **Mathematical Intuition:**
                
                - The Legendre symbol is used to determine whether a number has a square root modulo p (for an odd prime p).
                - It returns 1 if a is a quadratic residue modulo p, -1 if it is a non-residue, and 0 if a is divisible by p.
                
                **Preconditions:**
                
                - p must be a prime number.
                
                :param a: An integer a in ℤ.
                :param p: A Z object representing an odd prime modulus.
                :return: 1 if a is a quadratic residue modulo p, -1 if it is a non-residue, 0 if a is divisible by p.
                :raises Exception: If p is not an odd prime.
                """
                if int(cls.modulus) == 2:
                    raise Exception("legendre symbol not defined for p=2")
                symbol = pow(int(z.z), (int(cls.modulus) - 1) // 2, int(cls.modulus))
                return Z(-1) if symbol == int(cls.modulus) - 1 else Z(symbol)

            def __truediv__(self, b: Self) -> Self:
                if not b in self.__class__.units():
                    raise Exception("unsupported operation")
                return self * b.inverse()

    return Z_