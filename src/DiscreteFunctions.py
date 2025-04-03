import typing
import random
import math

class PrimalityTesting:
    r"""
    A collection of methods to test the primality of positive integers.
    
    **Preconditions:**
    
    - All methods assume that input integers are positive (n > 0).
    - For methods that use random sampling (Fermat and Miller–Rabin tests), it is assumed that n is sufficiently large (typically n > 3) to ensure a meaningful random base can be selected.
    """

    @staticmethod
    def coprime(a: int, b: int) -> bool:
        r"""
        Determine whether two integers are coprime, i.e., their greatest common divisor is 1.
        
        **Mathematical Intuition:**
        
        The integers a and b are coprime if:
        
        .. math::
            \gcd(a, b) = 1.
        
        This function verifies that no prime factor (other than 1) divides both a and b.
        
        **Preconditions:**
        
        - a and b are positive integers.
        
        :param a: A positive integer.
        :param b: A positive integer.
        :return: True if a and b are coprime; False otherwise.
        """
        return math.gcd(a, b) == 1

    @staticmethod
    def extended_euclidean(a: int, b: int) -> tuple[int, int]:
        r"""
        Compute the coefficients (x, y) that satisfy Bezout's identity:
        
        .. math::
            a \cdot x + b \cdot y = \gcd(a, b).
        
        **Mathematical Intuition:**
        
        This function returns integers x and y such that:
        
        .. math::
            a x + b y = \gcd(a, b),
        
        which is fundamental in solving Diophantine equations and computing modular inverses.
        
        **Preconditions:**
        
        - a and b are positive integers, not both zero.
        
        :param a: A positive integer.
        :param b: A positive integer.
        :return: A tuple (x, y) satisfying the equation above.
        """
        x1, y1 = 1, 0
        x2, y2 = 0, 1
    
        while b != 0:
            q = a // b
            a, b = b, a % b
            x1, x2 = x2, x1 - q * x2
            y1, y2 = y2, y1 - q * y2
    
        return (x1, y1)

    @staticmethod
    def eratosthenes_primality(n : int) -> bool:
        r"""
        Test if a number n is prime using a simple variant of the Sieve of Eratosthenes.
        
        **Mathematical Intuition:**
        
        An integer n > 1 is prime if it is not divisible by any integer in the range:
        
        .. math::
            2 \leq v \leq \sqrt{n}.
        
        This method checks divisibility only up to the integer square root of n.
        
        **Preconditions:**
        
        - n is an integer greater than 1.
        
        :param n: An integer > 1.
        :return: True if n is prime; False otherwise.
        """
        if n % 2 == 0:
            return n == 2
        for v in range(3, math.isqrt(n)+1, 2):
            if n % v == 0:
                return False
        return True

    @staticmethod
    def fermat_primality(n: int) -> bool:
        r"""
        Probabilistically test if a number is prime using Fermat's little theorem.
        
        **Mathematical Intuition:**
        
        Fermat's little theorem states that if n is prime and a is an integer with 1 < a < n, then:
        
        .. math::
            a^{n-1} \equiv 1 \pmod{n}.
        
        If this congruence fails for some a, then n is composite.
        
        **Preconditions:**
        
        - n is an integer greater than 3 (with 2 and 3 handled explicitly).
        - n is odd (as even numbers > 2 are composite).
        
        :param n: An integer greater than 3.
        :return: True if n passes Fermat's test (and is likely prime); False if n fails the test.
        """
        if n == 2 or n == 3:
            return True
        a = random.randint(2, math.isqrt(n))
        if not PrimalityTesting.coprime(a, n):
            return False
        return pow(a, n - 1, n) == 1

    @staticmethod
    def miller_rabin_primality(n: int) -> bool:
        r"""
        Test if a number is prime using the Miller–Rabin probabilistic test.
        
        **Mathematical Intuition:**
        
        Write n - 1 as:
        
        .. math::
            n - 1 = 2^s \cdot d,
        
        with d odd. Then for a random base a (2 ≤ a ≤ n − 2), compute:
        
        .. math::
            x = a^d \mod n.
        
        n is composite if:
        
        .. math::
            x \not\in \{1, n-1\} \quad \text{and} \quad x^{2^r} \not\equiv n-1 \, (\forall\, 0 \leq r < s).
        
        **Preconditions:**
        
        - n is an odd integer greater than 3 (with 2 and 3 handled explicitly).
        
        :param n: An odd integer > 3.
        :return: True if n passes the Miller–Rabin test (likely prime); False if n fails (composite).
        """
        if n == 2 or n == 3:
            return True
        s = int(math.log2((n-1) & -(n-1)))
        d = (n - 1) >> s
        a = random.randint(2, n-2)
        x = pow(a, d, n)
        if x != 1 and x != (n-1) and all([pow(a, (2 ** r)*d, n) != (n-1) for r in range(s)]):
            return False
        return True

class PrimeNumberTheorem:
    r"""
    Provides functions related to the prime number theorem (PNT) and approximations for the count of primes.
    
    **Preconditions:**
    
    - Input values are positive integers (typically x > 1).
    """

    @staticmethod
    def π(x: int) -> int:
        r"""
        Count the number of primes less than or equal to x.
        
        **Mathematical Intuition:**
        
        This function estimates:
        
        .. math::
            \pi(x) = \#\{ p \leq x \mid p \text{ is prime} \},
        
        by testing each integer n (2 ≤ n ≤ x) for primality.
        
        **Preconditions:**
        
        - x is an integer ≥ 2.
        
        :param x: A positive integer (x ≥ 2).
        :return: The count of primes ≤ x.
        """
        count = 0
        for n in range(2, x+1):
            if PrimalityTesting.fermat_primality(n):
                count += 1
        return count

    @staticmethod
    def prob(n: int) -> float:
        r"""
        Compute a heuristic probability that an integer n is prime.
        
        **Mathematical Intuition:**
        
        It uses a product over primes up to √n:
        
        .. math::
            P(n) \approx \prod_{f \leq \sqrt{n}} \frac{f-1}{f},
        
        which estimates the density of primes.
        
        **Preconditions:**
        
        - n is an integer ≥ 2.
        
        :param n: A positive integer.
        :return: A floating-point estimate of the probability that n is prime.
        """
        bound = math.floor(n ** 0.5)
        F = [((f-1)/f) for f in range(2, bound+1) if PrimalityTesting.fermat_primality(f)]
        return math.prod(F)

    @staticmethod
    def pnt_prime_count(x: int) -> float:
        r"""
        Estimate the number of primes up to x via summation.
        
        **Mathematical Intuition:**
        
        The prime number theorem suggests the prime density near n is ~1/ln(n), so:
        
        .. math::
            \pi(x) \approx \sum_{n=2}^{x} \frac{1}{\ln(n)}.
        
        **Preconditions:**
        
        - x is an integer > 1.
        
        :param x: A positive integer > 1.
        :return: A floating-point estimate of the number of primes ≤ x.
        """
        acc = 0
        for n in range(2, x+1):
            acc += 1 / math.log(n)
        return acc

    @staticmethod
    def approx_pnt_prime_count(x: int) -> float:
        r"""
        Approximate the number of primes up to x using the continuous form of the prime number theorem.
        
        **Mathematical Intuition:**
        
        The prime number theorem states:
        
        .. math::
            \pi(x) \sim \frac{x}{\ln(x)},
        
        so this function returns:
        
        .. math::
            \frac{x}{\ln(x)}.
        
        **Preconditions:**
        
        - x is an integer > 1.
        
        :param x: A positive integer > 1.
        :return: A floating-point approximation of the number of primes ≤ x.
        """
        return x / math.log(x)

class Factorization:
    r"""
    Provides methods for computing the divisors and prime factorization of a positive integer.
    
    **Preconditions:**
    
    - All functions assume that n is a positive integer (n ≥ 1).
    """
    @staticmethod
    def divisors(n: int) -> list[int]:
        r"""
        Compute all positive divisors of n.
        
        **Mathematical Intuition:**
        
        Divisors of n are numbers d such that:
        
        .. math::
            d \mid n.
        
        This function finds all such d by testing up to :math:`\sqrt{n}` and including their complementary divisors.
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        
        :param n: A positive integer.
        :return: A list of positive divisors of n.
        """
        res = [a for a in range(1, math.isqrt(n)+1) if n % a ==  0]
        return res + [n // f for f in res[::-1] if f * f != n]

    @staticmethod
    def prime_factors(n: int) -> list[int]:
        r"""
        Determine the distinct prime factors of n.
        
        **Mathematical Intuition:**
        
        A prime factor p of n satisfies:
        
        .. math::
            p \mid n \quad \text{and} \quad p \text{ is prime}.
        
        This function returns each prime factor once.
        
        **Preconditions:**
        
        - n is an integer greater than 1.
        
        :param n: A positive integer > 1.
        :return: A list of prime numbers dividing n.
        """
        return [a for a in Factorization.divisors(n)[1:] if PrimalityTesting.eratosthenes_primality(a)]

    @staticmethod
    def factorize(n: int) -> list[tuple[int, int]]:
        r"""
        Factorize n into its prime factors with exponents.
        
        **Mathematical Intuition:**
        
        For a positive integer n > 1, there exist primes :math:`p_1, p_2, \dots, p_k` and exponents :math:`e_1, e_2, \dots, e_k` such that:
        
        .. math::
            n = p_1^{e_1} \cdots p_k^{e_k}.
        
        This function returns a list of tuples :math:`(p_i, e_i)`.
        
        **Preconditions:**
        
        - n is an integer > 1.
        
        :param n: A positive integer > 1.
        :return: A list of tuples, where each tuple is (prime, exponent).
        """
        def helper(n,a):
            exp = 0
            while n % a == 0:
                n //= a
                exp += 1
            return exp
        return [(f, helper(n, f)) for f in Factorization.prime_factors(n)]

class ArithmeticFunctions:
    r"""
    Implements several arithmetic functions from number theory, including divisor functions and
    multiplicative functions such as Euler's totient function, Liouville's function, and the Möbius function.
    
    **Preconditions:**
    
    - The input n should be a positive integer (n ≥ 1).
    """

    @staticmethod
    def σ(n: int) -> int:
        r"""
        Compute the sum-of-divisors function σ(n).
        
        **Mathematical Intuition:**
        
        The function is defined by:
        
        .. math::
            \sigma(n) = \sum_{d \mid n} d,
        
        where the sum runs over all positive divisors d of n.
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        
        :param n: A positive integer.
        :return: The sum of all positive divisors of n.
        """
        return sum(Factorization.divisors(n))
    
    @staticmethod
    def σ_k(n: int, k: int) -> int:
        r"""
        Compute the generalized divisor function σₖ(n), which sums the k-th powers of all divisors of n.
        
        **Mathematical Intuition:**
        
        It is defined as:
        
        .. math::
            \sigma_k(n) = \sum_{d \mid n} d^k.
        
        For k = 1, this reduces to the standard sum-of-divisors function.
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        - k is an integer.
        
        :param n: A positive integer.
        :param k: An integer exponent.
        :return: The sum of each divisor of n raised to the power k.
        """
        return sum(d ** k for d in Factorization.divisors(n))

    @staticmethod
    def Τ(n: int) -> int:
        r"""
        Compute the divisor-counting function Τ(n), i.e., the total number of positive divisors of n.
        
        **Mathematical Intuition:**
        
        .. math::
            T(n) = \#\{d \mid n\}.
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        
        :param n: A positive integer.
        :return: The number of positive divisors of n.
        """
        return len(Factorization.divisors(n))

    @staticmethod
    def ω(n: int) -> int:
        r"""
        Compute the number of distinct prime factors of n.
        
        **Mathematical Intuition:**
        
        .. math::
            \omega(n) = \#\{p \mid n : p \text{ is prime}\}.
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        
        :param n: A positive integer.
        :return: The number of distinct prime factors of n.
        """
        return len(Factorization.prime_factors(n))

    @staticmethod
    def Ω(n: int) -> int:
        r"""
        Compute the total number of prime factors of n, counting multiplicities.
        
        **Mathematical Intuition:**
        
        If
        
        .. math::
            n = \prod_{i=1}^{k} p_i^{e_i},
        
        then
        
        .. math::
            \Omega(n) = \sum_{i=1}^{k} e_i.
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        
        :param n: A positive integer.
        :return: The total number of prime factors of n (with multiplicity).
        """
        return sum([a[1] for a in Factorization.factorize(n)])

    @staticmethod
    def ϕ(n: int) -> int:
        r"""
        Compute Euler's totient function ϕ(n), which counts the number of integers between 1 and n that are coprime to n.
        
        **Mathematical Intuition:**
        
        Euler's totient function is given by:
        
        .. math::
            \varphi(n) = n \prod_{p \mid n} \left(1 - \frac{1}{p}\right),
        
        where the product is over all distinct prime factors of n.
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        
        :param n: A positive integer.
        :return: The value of Euler's totient function ϕ(n).
        """
        return round(n * math.prod([1-(1/p) for p in Factorization.prime_factors(n)]))

    @staticmethod
    def λ(n: int) -> int:
        r"""
        Compute the Liouville function λ(n).
        
        **Mathematical Intuition:**
        
        The Liouville function is defined as:
        
        .. math::
            \lambda(n) = (-1)^{\Omega(n)},
        
        where Ω(n) is the total number of prime factors of n (with multiplicity).
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        
        :param n: A positive integer.
        :return: 1 if Ω(n) is even, -1 if Ω(n) is odd.
        """
        return (-1) ** ArithmeticFunctions.Ω(n)

    @staticmethod
    def μ(n: int) -> int:
        r"""
        Compute the Möbius function μ(n).
        
        **Mathematical Intuition:**
        
        The Möbius function is defined by:
        
        .. math::
            \mu(n) =
            \begin{cases}
            (-1)^{\omega(n)} & \text{if } n \text{ is square-free,} \\
            0 & \text{otherwise.}
            \end{cases}
        
        **Preconditions:**
        
        - n is an integer with n ≥ 1.
        
        :param n: A positive integer.
        :return: The Möbius function value: -1, 0, or 1.
        """
        if ArithmeticFunctions.ω(n) == ArithmeticFunctions.Ω(n):
            return (-1) ** ArithmeticFunctions.Ω(n)
        else:
            return 0

    @staticmethod
    def I(n: int) -> int:
        r"""
        The indicator function I(n), defined as:
        
        .. math::
            I(n) =
            \begin{cases}
            1 & \text{if } n = 1, \\
            0 & \text{otherwise.}
            \end{cases}
        
        **Preconditions:**
        
        - n is an integer.
        
        :param n: An integer.
        :return: 1 if n == 1, 0 otherwise.
        """
        if n == 1:
            return 1
        else:
            return 0

    @staticmethod
    def N(n: int) -> int:
        r"""
        The identity function N(n), which returns n itself.
        
        **Mathematical Intuition:**
        
        .. math::
            N(n) = n.
        
        This function acts as the natural number identity in various arithmetic convolutions.
        
        **Preconditions:**
        
        - n is an integer.
        
        :param n: An integer.
        :return: n.
        """
        return n
    
    @staticmethod
    def is_perfect(n: int) -> bool:
        r"""
        Determine whether n is a perfect number.
        
        **Mathematical Intuition:**
        
        A number n is perfect if:
        
        .. math::
            \sigma(n) = 2n,
        
        where σ(n) is the sum-of-divisors function.
        
        **Preconditions:**
        
        - n is a positive integer.
        
        :param n: A positive integer.
        :return: True if n is perfect; False otherwise.
        """
        return ArithmeticFunctions.σ(n) == 2 * n

    @staticmethod
    def is_square_free(n: int) -> bool:
        r"""
        Check if n is square-free (i.e., no prime factor appears with exponent greater than 1).
        
        **Mathematical Intuition:**
        
        n is square-free if:
        
        .. math::
            p^2 \nmid n \quad \text{for every prime } p.
        
        **Preconditions:**
        
        - n is a positive integer.
        
        :param n: A positive integer.
        :return: True if n is square-free; False otherwise.
        """
        if n < 1:
            return False
        if n == 1:
            return True

        for p in range(2, math.isqrt(n)+1):
            if PrimalityTesting.eratosthenes_primality(p) and n % (p*p) == 0:
                return False
        return True  