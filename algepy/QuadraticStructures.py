import typing
from typing import Callable, Self
import math
import sys
from .SingletonStructures import Z, R, Z_n, Z_mod_
from .DiscreteFunctions import ArithmeticFunctions, Factorization

class C:
    r"""
    A lightweight wrapper for Python's built-in complex numbers, used primarily for 
    embedding purposes in normalization routines for imaginary quadratic integer rings.
    
    **Mathematical Intuition:**
    
    Elements are represented as:
    
    .. math::
        c = a + bi, \quad a, b \in \mathbb{R},
    
    and arithmetic is performed using the standard operations on complex numbers.
    """
    def __init__(self, a, b) -> Self:
        try:
            a = float(a)
        except:
            try:
                a = int(a)
            except:
                raise Exception("a ¬∈ R")
        try:
            b = float(b)
        except:
            try:
                b = int(b)
            except:
                raise Exception("b ¬∈ R")
        self._value = complex(a, b)
    
    def __add__(self, c: Self) -> Self:
        res = self._value + c._value
        return C(res.real, res.imag)
    
    def __mul__(self, c) -> Self:
        if isinstance(c, C):
            res = self._value * c._value
            return C(res.real, res.imag)
        elif isinstance(c, R):
            res = self._value * c
            return C(res.real, res.imag)
        return NotImplemented
    
    def __rmul__(self, c) -> Self:
        return self.__mul__(c)
    
    def __neg__(self) -> Self:
        res = -self._value
        return C(res.real, res.imag)
    
    def __sub__(self, c: Self) -> Self:
        res = self._value - c._value
        return C(res.real, res.imag)
    
    def conjugate(self) -> Self:
        r"""
        Returns the complex conjugate of the element.
        
        **Mathematical Intuition:**
        
        For an element :math:`c = a + bi`, the conjugate is defined as:
        
        .. math::
            \overline{c} = a - bi.
        
        This operation is used in the computation of norms and is fundamental in 
        normalizing elements of imaginary quadratic rings.
        """
        res = self._value.conjugate()
        return C(res.real, res.imag)
    
    def __truediv__(self, c: Self) -> Self:
        res = self._value / c._value
        return C(res.real, res.imag)
    
    def __eq__(self, c: Self) -> bool:
        return self._value == c._value
    
    def __pow__(self, e: Z) -> Self:
        res = self._value ** int(e)
        return C(res.real, res.imag)
    
    def __repr__(self) -> str:
        return str(self._value.real) + " + " + str(self._value.imag) + "i"
    
    def to_tuple(self) -> tuple[float, float]:
        return (self._value.real, self._value.imag)
    
    @staticmethod
    def zero() -> Self:
        r"""
        Returns the additive identity, :math:`0 + 0i`, in the complex wrapper.

        **Mathematical Intuition:**
        
        .. math::
            0 \in \mathbb{C}, \quad \text{such that } z + 0 = z.
        
        :return: The zero element as a C object.
        """
        return C(0, 0)
    
    @staticmethod
    def mul_id() -> Self:
        r"""
        Returns the multiplicative identity, :math:`1 + 0i`, in the complex wrapper.
        
        **Mathematical Intuition:**
        
        .. math::
            1 \in \mathbb{C}, \quad \text{such that } z \cdot 1 = z.
        
        :return: The identity element as a C object.
        """
        return C(1, 0)
    
    @staticmethod
    def norm(c: Self) -> R:
        r"""
        Computes the norm of a complex number.
        
        **Mathematical Intuition:**
        
        The norm is defined as:
        
        .. math::
            \|c\| = \sqrt{a^2 + b^2} \quad \text{for } c = a+bi.
        
        :param c: A complex element.
        :return: The norm as an element of R.
        """
        return R(((c._value.real ** 2) + (c._value.imag ** 2)) ** 0.5)

class Q:
    r"""
    Represents an element of the quadratic rational field :math:`\mathbb{Q}(\sqrt{d})` in standard form.
    
    An element is stored in reduced form as:
    
    .. math::
        q = \frac{n}{d}, \quad \text{with } \gcd(n,d)=1 \text{ and } d > 0.
    
    **Mathematical Intuition:**
    
    This representation ensures each element is uniquely determined by its numerator and denominator,
    analogous to a rational number.
    """
    def __init__(self, n, d) -> Self:
        
        try:
            if int(d) == 0:
                raise Exception()
            
            gcd = Z(math.gcd(int(n), int(d)))
            self.n = Z.abs((n if isinstance(n, Z) else Z(n)) // gcd) * (Z(-1) if int(n) * int(d) < 0 else Z.mul_id())
            self.d = Z.abs((d if isinstance(d, Z) else Z(d)) // gcd)
        except:
            raise Exception("n ¬∈ Z ∨ d ¬∈ Z ∨ d = 0")

    def __add__(self, q: Self) -> Self:
        lcm = Z(math.lcm(int(self.d), int(q.d)))
        return Q(self.n * (lcm // self.d) + q.n * (lcm // q.d), lcm)

    def __mul__(self, q: Self) -> Self:
        if isinstance(q, Q):
            return Q(self.n * q.n, self.d * q.d)
        elif isinstance(q, Z):
            return Q(self.n * q, self.d * q)
        return NotImplemented
    
    def __rmul__(self, q) -> Self:
        return self.__mul__(q)

    def __neg__(self) -> Self:
        return Q(-self.n, self.d)

    def __sub__(self, q: Self) -> Self:
        return self + -q

    def __truediv__(self, q: Self) -> Self:
        return self * Q(q.d, q.n)

    def __eq__(self, q: Self) -> bool:
        return self.n == q.n and self.d == q.d

    def __repr__(self) -> str:
        return (str(self.n) + "/" + str(self.d))

    def __int__(self) -> int:
        if self.d != Z(1):
            raise Exception("not an integer")
        return int(self.n)

    def __float__(self) -> float:
        return int(self.n) / int(self.d)
    
    def reciprocal(self) -> Self:
        r"""
        Returns the multiplicative inverse of this rational element.
        
        **Mathematical Intuition:**
        
        For a nonzero :math:`\frac{n}{d}`, its reciprocal is:
        
        .. math::
            \left(\frac{n}{d}\right)^{-1} = \frac{d}{n}.
        
        **Preconditions:**
        
        - :math:`\frac{n}{d} \neq 0`.
        """
        return Q(self.d, self.n)
    
    def to_tuple(self) -> tuple[int, int]:
        return (self.n, self.d)
    
    @staticmethod
    def zero() -> Self:
        r"""
        Returns the additive identity :math:`\frac{0}{1}` in this rational form.
        
        **Mathematical Intuition:**
        
        .. math::
            0 \in \mathbb{Q}, \quad \text{such that } q + 0 = q.
        
        :return: The zero element as a Q object.
        """
        return Q(0, 1)
    
    @staticmethod
    def mul_id() -> Self:
        r"""
        Returns the multiplicative identity :math:`\frac{1}{1}`.
        
        **Mathematical Intuition:**
        
        .. math::
            1 \in \mathbb{Q}, \quad \text{such that } q \cdot 1 = q.
        
        :return: The identity element as a Q object.
        """
        return Q(1, 1)

    @staticmethod
    def abs(q: Self) -> Self:
        return Q(Z.abs(q.n), Z.abs(q.d))

class QuadInt:
    r"""
    Represents an element of a quadratic integer ring (e.g. :math:`\mathbb{Z}[\sqrt{d}]` or its half-integral variant).
    
    An element is represented as:
    
    .. math::
        a + b\omega,
    
    where :math:`\omega = \sqrt{d}` (or :math:`\omega = \frac{1+\sqrt{d}}{2}` when :math:`d \equiv 1 \mod 4`).
    
    **Mathematical Intuition:**
    
    Quadratic integers form an integral domain, and when the ring is a UFD, every element factors uniquely 
    (up to units). The operations (addition, multiplication, etc.) follow the usual algebraic rules.
    """
    def __init__(self, a, b, d):
        try:
            self.gen = d if isinstance(d, Z) else Z(d)
            self.a = a if isinstance(a, Z) else Z(a)
            self.b = b if isinstance(b, Z) else Z(b)
        except:
            raise Exception("a ¬∈ Z ∨ b ¬∈ Z ∨ d invalid")
        
    def __eq__(self, z: Self) -> bool:
        return self.a == z.a and self.b == z.b and self.gen == z.gen
    
    def __add__(self, z: Self) -> Self:
        if not self.gen == z.gen:
            raise Exception("unsupported operation")
        return self.__class__(self.a + z.a, self.b + z.b, self.gen)
    
    def __neg__(self) -> Self:
        return self.__class__(-self.a, -self.b, self.gen)
    
    def __sub__(self, z: Self) -> Self:
        if not self.gen == z.gen:
            raise Exception("unsupported operation")
        return self.__class__(self.a + -z.a, self.b + -z.b, self.gen)
        
    def __mul__(self, z) -> Self:
        if isinstance(z, QuadInt) and self.gen == z.gen:
            if self.gen % Z(4) == Z(1):
                return self.__class__(self.a * z.a + (self.b * z.b * (self.gen // Z(4))), self.a * z.b + self.b * z.a + self.b * z.b, self.gen)
            else:
                return self.__class__(self.a * z.a + (self.b * z.b * self.gen), self.a * z.b + self.b * z.a, self.gen)
        elif isinstance(z, Z):
            return self.__class__(self.a * z, self.b * z, self.gen)
        return NotImplemented
    
    def __rmul__(self, z) -> Self:
        return self.__mul__(z)
    
    def conjugate(self) -> Self:
        r"""
        Returns the conjugate of the quadratic integer.
        
        **Mathematical Intuition:**
        
        For an element

        .. math::
            a + b\sqrt{d},

        the conjugate is defined as:

        .. math::
            a - b\sqrt{d}.
        
        In the half-integral case (when :math:`d \equiv 1 \mod 4`), the conjugation is adjusted accordingly.
        """
        return self.__class__(self.a + self.b, -self.b, self.gen) if self.gen % Z(4) == Z(1) else self.__class__(self.a, -self.b, self.gen)
    
    @staticmethod
    def norm(z: Self) -> Z:
        r"""
        Computes the norm of a quadratic integer.
        
        **Mathematical Intuition:**
        
        For :math:`z = a + b\sqrt{d}` (or its half-integral variant), the norm is given by:
        
        .. math::
            N(z) = z \cdot \overline{z} = a^2 - d\,b^2.
        
        The norm is an integer and is multiplicative.
        """
        product = z * z.conjugate()
        return product.a
    
    @staticmethod
    def trace(z: Self) -> Z:
        r"""
        Computes the trace of a quadratic integer.
        
        **Mathematical Intuition:**
        
        The trace is defined as:
        
        .. math::
            \operatorname{Tr}(z) = z + \overline{z} = 2a \quad \text{(for standard form)},
        
        and is an integer.
        """
        sum = z + z.conjugate()
        return sum.a
    
    def inverse(self) -> Self:
        r"""
        Computes the inverse of a quadratic integer, when it exists.
        
        **Mathematical Intuition:**
        
        For a unit :math:`z` (i.e., one with norm ±1), the inverse is given by:
        
        .. math::
            z^{-1} = \frac{\overline{z}}{N(z)}.
        
        **Preconditions:**
        
        - :math:`z` must be a unit (i.e., :math:`|N(z)| = 1`).
        
        :return: The inverse of :math:`z`.
        :raises Exception: If :math:`z` is not a unit.
        """
        norm = QuadInt.norm(self)
        if int(norm) not in {-1,1}:
            raise Exception("unsupported operation")
        return self.conjugate() * norm

    def __pow__(self, e: Z) -> Self:
        r"""
        Exponentiates the quadratic integer using exponentiation by squaring.
        
        **Mathematical Intuition:**
        
        For any integer :math:`e`, the power :math:`z^e` is computed in :math:`O(\log |e|)` multiplications.
        
        **Preconditions:**
        
        - :math:`e` is a nonnegative integer.
        """
        if e == Z.zero():
            return self.__class__(1, 0, self.gen)
        
        result = self.__class__(1, 0, self.gen)
        base = self
        while not e == Z.zero():
            if e % Z(2) == Z(1):
                result = result * base
            base = base * base
            e = e // Z(2)

        return result

    def __repr__(self) -> str:
        omega = "((1+√"+str(self.gen)+")/2)" if self.gen % Z(4) == Z(1) else "(√"+str(self.gen)+")"
        return str(self.a) + " + " + str(self.b) + omega
    
def QuadIntRing(d: int, force_ufd: bool = False) -> Callable[[int, int], QuadInt]:
    r"""
    Constructs the quadratic integer ring with generator :math:`\omega` determined by d.
    
    **Mathematical Intuition:**
    
    - For a square-free integer :math:`d \neq 0,1`, the quadratic integer ring is given by:
      
      .. math::
          \mathbb{Z}[\sqrt{d}] \quad \text{if } d \not\equiv 1 \pmod{4},
      
      or
      
      .. math::
          \mathbb{Z}\left[\frac{1+\sqrt{d}}{2}\right] \quad \text{if } d \equiv 1 \pmod{4}.
    
    - Additional operations (such as division in the Euclidean case, normalization, and factorization) 
      are provided when the ring is known to be norm-Euclidean or a UFD.
    
    **Preconditions:**
    
    - :math:`d` must be a square-free integer with :math:`d \not\in \{0,1\}`.
    
    :param d: A square-free integer (not 0 or 1) representing the parameter of the quadratic field.
    :param known_ufd: Boolean flag indicating if the ring is known to be a UFD.
    :return: A callable class for constructing elements of the quadratic integer ring.
    :raises Exception: If d is invalid.
    """
    if int(d) in {0,1} or not ArithmeticFunctions.is_square_free(abs(d)):
                raise Exception("d invalid")
    
    class Quad_Z_(QuadInt):
        """
        NOTE: This is the subclass constructed by function QuadIntRing().
        """
        omega = Z(d)
        is_euclidean_domain = False
        is_ufd = force_ufd

        def __init__(self, a, b, d=None):
            super().__init__(a, b, self.omega)

        @classmethod
        def zero(cls) -> Self:
            r"""
            Returns the additive identity :math:`0` in the quadratic integer ring.
            
            **Mathematical Intuition:**
            
            .. math::
                0 + 0\omega,
            
            is the unique element satisfying :math:`z + 0 = z`.
            """
            return cls(0,0)
        
        @classmethod
        def mul_id(cls) -> Self:
            r"""
            Returns the multiplicative identity :math:`1` in the quadratic integer ring.
            
            **Mathematical Intuition:**
            
            .. math::
                1 + 0\omega,
            
            is the unique element satisfying :math:`z \cdot 1 = z`.
            """
            return cls(1,0)
        
        if d in {-11, -7, -3, -2, -1, 2, 3, 5, 6, 7, 11, 13, 17, 19, 21, 29, 33, 37, 41, 57, 73}:
            is_euclidean_domain = True
            is_ufd = True

            @classmethod
            def internal_div(cls, self: Self, z: Self) -> Self:
                r"""
                Computes the floor division of two quadratic integers.
                
                **Mathematical Intuition:**
                
                For norm-Euclidean domains, given :math:`w` and a nonzero :math:`z`, there exists a quotient :math:`q`
                and remainder :math:`r` such that:
                
                .. math::
                    w = zq + r \quad \text{with} \quad N(r) < N(z).
                
                This method returns the quotient :math:`q`.
                
                **Preconditions:**
                
                - :math:`z` must be nonzero and belong to the same ring.
                """
                if z.a == Z.zero() and z.b == Z.zero():
                    raise ZeroDivisionError("Division by zero")
                if not self.gen == z.gen:
                    raise Exception("unsupported operation")

                norm_z = QuadInt.norm(z)
                numerator = self * z.conjugate()

                x = int(numerator.a) / int(norm_z)
                y = int(numerator.b) / int(norm_z)
                d_val = int(self.gen)
                if d_val % 4 != 1:

                    q1 = round(x)
                    q2 = round(y)
                    return cls(q1, q2)
                else:
                    b_candidate = round(2 * y)
                    a_candidate = round(x - (b_candidate / 2))
                    return cls(a_candidate, b_candidate)
            
            @classmethod
            def internal_mod(cls, self: Self, z: Self) -> Self:
                return self - (z * (cls.internal_div(self, z)))
            
            def __floordiv(self, z: Self) -> Self:
                return self.__class__.internal_div(self, z)
            
            def __mod__(self, z: Self) -> Self:
                return self.__class__.internal_mod(self, z)
            
            @classmethod
            def gcd(cls, w: Self, z: Self) -> Self:
                r"""
                Computes the greatest common divisor (gcd) of two quadratic integers via the Euclidean algorithm.
                
                **Mathematical Intuition:**
                
                In a Euclidean domain, for :math:`w` and :math:`z`, the gcd is computed as:
                
                .. math::
                    \gcd(w, z) = \gcd(z, w \mod z),
                
                iterated until the remainder is zero.
                
                **Preconditions:**
                
                - :math:`w` and :math:`z` belong to the same ring.
                """
                # w = w.normalize()
                # z = z.normalize()
                while not z == cls.zero():
                    r = cls.internal_mod(w, z)
                    w = z
                    z = r
                return w

        if is_euclidean_domain or is_ufd:
            is_ufd = True
            
            if d > 1:
                @classmethod
                def fundamental_unit(cls) -> Self:
                    r"""
                    Computes the fundamental unit of the real quadratic ring.
                    
                    **Mathematical Intuition:**
                    
                    For :math:`\mathbb{Z}[\sqrt{d}]` (or its half-integral variant), the fundamental unit 
                    :math:`\varepsilon` is the smallest unit :math:`> 1` (in absolute value) satisfying:
                    
                    .. math::
                        \varepsilon = x + y\sqrt{d} \quad \text{with} \quad x^2 - d\,y^2 = \pm 1.
                    
                    The method uses a continued fraction approach.
                    
                    **Preconditions:**
                    
                    - :math:`d > 1` and square-free.
                    """
                    if int(cls.omega) % 4 == 1:
                        return cls(0, 1)
                    a0 = int(math.floor(math.sqrt(d)))
                    m = 0
                    d0 = 1
                    a = a0
                    period = []
                    while True:
                        m = a * d0 - m
                        d0 = (int(cls.omega) - m * m) // d0
                        a = (a0 + m) // d0
                        period.append(a)
                        if a == 2 * a0:
                            break
                    k = len(period)
                    target = 1 if (k % 2 == 1) else k
                    p = [1, a0]
                    q = [0, 1]
                    for i in range(1, target):
                        a_i = period[(i - 1) % k]
                        p_next = a_i * p[-1] + p[-2]
                        q_next = a_i * q[-1] + q[-2]
                        p.append(p_next)
                        q.append(q_next)
                    x, y = p[target], q[target]
                    return cls(x, y)
                
                epsilon = None
                
                @classmethod
                def get_fundamental_unit(cls) -> Self:
                    r"""
                    Returns the fundamental unit of the ring.
                    
                    **Mathematical Intuition:**
                    
                    The fundamental unit is the generator of the infinite part of the unit group in a real quadratic field.
                    """
                    if cls.epsilon is None:
                        cls.epsilon = cls.fundamental_unit()
                    return cls.epsilon
                
                @classmethod
                def embed(cls, z: Self) -> R:
                    r"""
                    Embeds a quadratic integer into the real numbers.
                    
                    **Mathematical Intuition:**
                    
                    For :math:`z = a + b\sqrt{d}`, the embedding is:
                    
                    .. math::
                        \sigma(z) = a + b\sqrt{d} \in \mathbb{R}.
                    """
                    if int(cls.omega) % 4 == 1:
                        return R(int(z.a) + (int(z.b) / 2.0) + ((int(z.b) / 2.0) * math.sqrt(int(cls.omega))))
                    else:
                        return R(int(z.a) + (int(z.b) * math.sqrt(int(cls.omega))))
                    
                @classmethod
                def unit_power(cls, u: Self, n: int) -> Self:
                    if n >= 0:
                        return u.__pow__(Z(n))
                    else:
                        return u.__pow__(-Z(n)).inverse()
            
                def normalize(self) -> Self:
                    r"""
                    Returns the canonical representative (via normalization) of the associate class.
                    
                    **Mathematical Intuition:**
                    
                    For real quadratic rings, every element :math:`z` has associates of the form 
                    :math:`\varepsilon^n z` (with :math:`\varepsilon` the fundamental unit). Normalization 
                    selects the unique representative :math:`z_{\mathrm{can}}` such that:
                    
                    .. math::
                        1 \leq |\sigma(z_{\mathrm{can}})| < |\sigma(\varepsilon)|.
                    
                    **Preconditions:**
                    
                    - :math:`d > 1`.
                    """
                    fundamental = self.__class__.get_fundamental_unit()
                    sigma_epsilon = abs(float(self.__class__.embed(fundamental)))
                    
                    sigma_z = abs(float(self.__class__.embed(self)))
                    
                    n = math.floor(math.log(sigma_z) / math.log(sigma_epsilon))
                    
                    candidate = self * self.__class__.unit_power(fundamental, -n)
                    emb = abs(float(self.__class__.embed(candidate)))
                    
                    while emb < 1:
                        candidate = candidate * fundamental
                        emb = abs(float(self.__class__.embed(candidate)))
                    
                    while emb >= sigma_epsilon:
                        candidate = candidate * fundamental.inverse()
                        emb = abs(float(self.__class__.embed(candidate)))
                    
                    return candidate
                
            else:
                @classmethod
                def imaginary_units(cls) -> list[Self]:
                    r"""
                    Returns the finite set of units for an imaginary quadratic ring.
                    
                    **Mathematical Intuition:**
                    
                    For an imaginary quadratic ring, the unit group is finite. For example:
                    
                    - If :math:`d = -1`, the units are :math:`\{1, -1, i, -i\}`.
                    - If :math:`d = -3`, the units are the 6th roots of unity.
                    
                    :param d: The parameter of the quadratic field.
                    :return: A set of quadratic integer units.
                    """
                    if int(cls.omega) == -1:
                        return [
                            cls(1, 0),
                            cls(-1, 0),
                            cls(0, 1),
                            cls(0, -1)
                        ]
                    elif int(cls.omega) == -3:
                        omega = cls(0, 1)
                        omega2 = omega - cls(1, 0)
                        return [
                            cls(1, 0), 
                            cls(-1, 0), 
                            omega, 
                            -omega, 
                            omega2, 
                            -omega2]
                    else:
                        return [
                            cls(1, 0),
                            cls(-1, 0)
                        ]

                units = None

                @classmethod
                def get_imaginary_units(cls) -> list[Self]:
                    r"""
                    Returns the set of units for the imaginary quadratic ring.
                    """
                    if cls.units is None:
                        cls.units = cls.imaginary_units()
                    return cls.units
                
                @classmethod
                def embed(cls, z: Self) -> C:
                    r"""
                    Embeds a quadratic integer from an imaginary quadratic ring into :math:`\mathbb{C}`.
                    
                    **Mathematical Intuition:**
                    
                    For :math:`z = a + b\sqrt{d}` with :math:`d < 0`, the natural embedding is:
                    
                    .. math::
                        \sigma(z) = a + i\,b\sqrt{-d} \in \mathbb{C}.
                    """
                    if int(cls.omega) % 4 == 1:
                        return C(int(z.a) + int(z.b) / 2.0, (int(z.b) / 2.0) * math.sqrt(-int(cls.omega)))
                    else:
                        return C(int(z.a), int(z.b) * math.sqrt(-int(cls.omega)))
                
                def normalize(self) -> Self:
                    r"""
                    Returns the canonical representative (via normalization) of the associate class in an imaginary quadratic ring.
                    
                    **Mathematical Intuition:**
                    
                    Since the unit group is finite, normalization is achieved by selecting the associate with 
                    minimal lexicographic order of its complex embedding.
                    
                    **Preconditions:**
                    
                    - :math:`d < 1` (i.e., an imaginary quadratic ring).
                    """
                    d = int(self.gen)
                    
                    U = self.__class__.get_imaginary_units()
                    
                    best_candidate = None
                    best_value = None
                    
                    for u in U:
                        candidate = u * self
                        candidate_val = self.__class__.embed(candidate).to_tuple()
                        
                        if best_candidate is None or candidate_val < best_value:
                            best_candidate = candidate
                            best_value = candidate_val
                            
                    return best_candidate
                
            @classmethod
            def factorize(cls, z: Self) -> list[tuple[Self, Z]]:
                r"""
                Factorizes a quadratic integer into irreducible factors (up to units).
                
                **Mathematical Intuition:**
                
                In a UFD, every nonzero non-unit element can be written uniquely (up to units) as:
                
                .. math::
                    z = \pi_1^{e_1} \pi_2^{e_2} \cdots \pi_k^{e_k}.
                
                This method computes such a factorization by using the norm to guide a trial division approach.
                It terminates because the norm is a positive integer that strictly decreases with each factor extraction.
                
                **Preconditions:**
                
                - :math:`z` must be a nonzero element of the quadratic integer ring.

                :return: A list of tuples, each consisting of an irreducible factor and its exponent.
                """
                def is_unit(x: Self) -> bool:
                    return abs(int(QuadInt.norm(x))) == 1

                normalized_z = z.normalize()
                current = normalized_z
                factors: list[tuple[Self, Z]] = []
                if is_unit(current):
                    return factors

                n_int = abs(int(QuadInt.norm(current)))
                prime_factors = Factorization.prime_factors(n_int)

                for p in prime_factors:
                    bound = int(math.sqrt(p)) + 2
                    factor_found = True
                    while factor_found and not is_unit(current):
                        factor_found = False
                        found = False
                        for a in range(-bound, bound + 1):
                            for b in range(-bound, bound + 1):
                                if a == 0 and b == 0:
                                    continue
                                candidate = cls(a, b)
                                cand_norm = abs(int(QuadInt.norm(candidate)))
                                if cand_norm not in {p, p * p}:
                                    continue
                                if is_unit(candidate):
                                    continue
                                for cand in [candidate, -candidate]:
                                    try:
                                        if current.__class__.internal_mod(current, cand) == cls.zero():
                                            exp = Z(0)
                                            while current.__class__.internal_mod(current, cand) == cls.zero():
                                                current = current.__class__.internal_div(current, cand)
                                                exp = exp + Z(1)
                                                current = current.normalize()
                                            factors.append((cand.normalize(), exp))
                                            factor_found = True
                                            found = True
                                            break
                                    except Exception:
                                        continue
                                if found:
                                    break
                            if found:
                                break
                        if not factor_found and not is_unit(current):
                            bound += 1

                if not is_unit(current):
                    factors.append((current.normalize(), Z(1)))

                prod = cls.mul_id()
                for fac, exp in factors:
                    for _ in range(int(exp)):
                        prod = prod * fac
                prod = prod.normalize()

                if prod != normalized_z:
                    u = normalized_z.__class__.internal_div(normalized_z, prod)
                    if factors:
                        f, e = factors[0]
                        if int(e) > 1:
                            factors[0] = (f * u, Z(1))
                            if int(e) - 1 > 0:
                                factors.insert(1, (f, e - Z(1)))
                        else:
                            factors[0] = (f * u, Z(1))
                return factors
    Quad_Z_.__module__ = __name__
    if 'sphinx' in sys.modules:
        globals()[f'Quad_Z_'] = Quad_Z_

    return Quad_Z_
                    
Quad_Z_ = QuadIntRing(5)

class QuadRat:
    r"""
    Represents an element of a quadratic rational field :math:`\mathbb{Q}(\sqrt{d})` in standard form.
    
    Each element is represented as:
    
    .. math::
        a + b\sqrt{d}, \quad a,b \in \mathbb{Q},
    
    where the representation is assumed to be in lowest terms.
    
    **Mathematical Intuition:**
    
    This class models the field :math:`\mathbb{Q}(\sqrt{d})`, in which every nonzero element is invertible.
    Arithmetic operations follow the usual rules for addition, multiplication, and inversion in a field.
    """
    def __init__(self, a, b, d):
        try:
            self.gen = d if isinstance(d, Z) else Z(d)
            if isinstance(a, Q) and isinstance(b, Q):
                self.a = a
                self.b = b
            else:
                raise Exception()
        except:
            raise Exception("a ¬∈ Q ∨ b ¬∈ Q ∨ d invalid")
        
    def __eq__(self, q: Self) -> bool:
        return self.a == q.a and self.b == q.b and self.gen == q.gen
    
    def __add__(self, q: Self) -> Self:
        if not self.gen == q.gen:
            raise Exception("unsupported operation")
        return self.__class__(self.a + q.a, self.b + q.b, self.gen)
    
    def __neg__(self) -> Self:
        return self.__class__(-self.a, -self.b, self.gen)
    
    def __sub__(self, q: Self) -> Self:
        if not self.gen == q.gen:
            raise Exception("unsupported operation")
        return self + -q
    
    def __mul__(self, q: Self) -> Self:
        if isinstance(q, QuadRat) and self.gen == q.gen:
            return self.__class__(self.a * q.a + (self.b * q.b * self.gen), self.a * q.b + self.b * q.a, self.gen)
        elif isinstance(q, Q):
            return self.__class__(self.a * q, self.b * q, self.gen)
        return NotImplemented
    
    def __rmul__(self, q) -> Self:
        return self.__mul__(q)
    
    def conjugate(self) -> Self:
        r"""
        Returns the conjugate of a quadratic rational element.
        
        **Mathematical Intuition:**
        
        For an element :math:`q = a + b\sqrt{d} \in \mathbb{Q}(\sqrt{d})`, the conjugate is:
        
        .. math::
            \overline{q} = a - b\sqrt{d}.
        
        This operation is essential for computing the norm and inversion.
        """
        return self.__class__(self.a, -self.b, self.gen)
    
    @staticmethod
    def norm(q: Self) -> Q:
        r"""
        Computes the norm of a quadratic rational element.
        
        **Mathematical Intuition:**
        
        For :math:`q = a + b\sqrt{d}`, the norm is given by:
        
        .. math::
            N(q) = (a+b\sqrt{d})(a-b\sqrt{d}) = a^2 - d\,b^2.
        
        The norm is a rational number.
        """
        product = q * q.conjugate()
        return product.a
    
    @staticmethod
    def trace(q: Self) -> Z:
        r"""
        Computes the trace of a quadratic rational element.
        
        **Mathematical Intuition:**
        
        The trace is defined by:
        
        .. math::
            \operatorname{Tr}(q) = q + \overline{q} = 2a.
        
        It is an element of ℤ.
        """
        sum = q + q.conjugate()
        return sum.a
    
    def inverse(self) -> Self:
        r"""
        Computes the inverse of a nonzero quadratic rational element.
        
        **Mathematical Intuition:**
        
        Given :math:`q = a + b\sqrt{d}` with nonzero norm, its inverse is:
        
        .. math::
            q^{-1} = \frac{\overline{q}}{N(q)} = \frac{a-b\sqrt{d}}{a^2-d\,b^2}.
        
        **Preconditions:**
        
        - :math:`q` must be nonzero.
        """
        norm = QuadRat.norm(self)
        if norm == Q.zero():
            raise Exception("unsupported operation")
        return norm.reciprocal() * self.conjugate()
    
    def __truediv__(self, q: Self) -> Self:
        return self * q.inverse()
    
    def __pow__(self, e: Z) -> Self:
        r"""
        Exponentiates the quadratic rational element using exponentiation by squaring.
        
        **Mathematical Intuition:**
        
        For any integer :math:`e`, :math:`q^e` is computed efficiently in :math:`O(\log |e|)` steps.
        
        **Preconditions:**
        
        - :math:`e` is an integer.
        """
        if e == Z.zero():
            return self.__class__(Q.mul_id(), Q.zero(), self.gen)
        
        result = self.__class__(Q.mul_id(), Q.zero(), self.gen)
        base = self
        while not e == Z.zero():
            if e % Z(2) == Z(1):
                result = result * base
            base = base * base
            e = e // Z(2)

        return result
    
    def __repr__(self) -> str:
        omega = "(√"+self.gen+")"
        return str(self.a) + " + " + str(self.b) + omega

def QuadRatField(d: int) -> Callable[[Q, Q], QuadRat]:
    r"""
    Constructs the quadratic rational field :math:`\mathbb{Q}(\sqrt{d})` in standard form.
    
    **Mathematical Intuition:**
    
    For a square-free integer :math:`d \neq 0,1`, the field is:
    
    .. math::
        \mathbb{Q}(\sqrt{d}) = \{\, a + b\sqrt{d} : a, b \in \mathbb{Q} \,\}.
    
    The field is represented in standard form; all arithmetic operations follow the usual field properties.
    
    **Preconditions:**
    
    - :math:`d` must be a nonzero square-free integer, and :math:`d \neq 1`.
    
    :param d: A square-free integer (not 0 or 1).
    :return: A callable class for constructing elements of :math:`\mathbb{Q}(\sqrt{d})`.
    :raises Exception: If d is invalid.
    """
    if int(d) in {0,1} or not ArithmeticFunctions.is_square_free(d):
                raise Exception("d invalid")
    
    class Quad_Q_(QuadRat):
        """
        NOTE: This is the subclass constructed by function QuadRatField().
        """
        omega = Z(d)

        def __init__(self, a, b, d=None):
            super().__init__(a, b, self.omega)

        @classmethod
        def integer_ring(cls):
            return QuadIntRing(int(cls.omega))
        
        @classmethod
        def discriminant(cls) -> Z:
            r"""
            Computes the field discriminant of :math:`\mathbb{Q}(\sqrt{d})`.
            
            **Mathematical Intuition:**
            
            The field discriminant is given by:
            
            .. math::
                \Delta =
                \begin{cases}
                d, & \text{if } d \equiv 1 \mod 4, \\
                4d, & \text{otherwise.}
                \end{cases}
            
            :return: The discriminant as a Z object.
            """
            return cls.omega if cls.omega % Z(4) == Z(1) else Z(4) * cls.omega
        
        if d > 1:
            @classmethod
            def embed(cls, q: Self) -> R:
                r"""
                Embeds an element :math:`q = a+b\sqrt{d}` into :math:`\mathbb{R}`.
                
                **Mathematical Intuition:**
                
                The standard embedding is:
                
                .. math::
                    \sigma(q) = a + b\sqrt{d}.
                
                :param q: An element of :math:`\mathbb{Q}(\sqrt{d})`.
                :return: The real number :math:`\sigma(q)`.
                """
                return R(float(q.a) + float(q.b) * math.sqrt(d))

            @classmethod
            def embed_prime(cls, q: Self) -> R:
                r"""
                Provides the second real embedding for :math:`q = a+b\sqrt{d}` into :math:`\mathbb{R}`.
                
                **Mathematical Intuition:**
                
                The second embedding is:
                
                .. math::
                    \sigma'(q) = a - b\sqrt{d}.
                
                :param q: An element of :math:`\mathbb{Q}(\sqrt{d})`.
                :return: The real number :math:`\sigma'(q)`.
                """
                return R(float(q.a) - float(q.b) * math.sqrt(d))
        else:
            @classmethod
            def embed(cls, q: Self) -> C:
                r"""
                Embeds an element :math:`q = a+b\sqrt{d}` (with :math:`d < 0`) into :math:`\mathbb{C}`.
                
                **Mathematical Intuition:**
                
                The standard complex embedding is:
                
                .. math::
                    \sigma(q) = a + i\,b\sqrt{-d}.
                
                :param q: An element of :math:`\mathbb{Q}(\sqrt{d})`.
                :return: A C object representing the complex number :math:`\sigma(q)`.
                """
                return C(float(q.a), float(q.b) * math.sqrt(-d))

            @classmethod
            def embed_prime(cls, q: Self) -> C:
                r"""
                Provides the second complex embedding for :math:`q = a+b\sqrt{d}` (with :math:` d < 0`).
                
                **Mathematical Intuition:**
                
                The second embedding is:
                
                .. math::
                    \sigma'(q) = a - i\,b\sqrt{-d}.
                
                :param q: An element of :math:`\mathbb{Q}(\sqrt{d})`.
                :return: A C object representing :math:`\sigma'(q)`.
                """
                return C(float(q.a), -float(q.b) * math.sqrt(-d))
            
    Quad_Q_.__module__ = __name__
    if 'sphinx' in sys.modules:
        globals()[f'Quad_Q_'] = Quad_Q_

    return Quad_Q_

Quad_Q_ = QuadRatField(5)
        

        

        
