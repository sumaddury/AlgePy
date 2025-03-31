import typing
import random
import math

class PrimalityTesting:
    """
    A collection of methods to test the primality of integers.
    
    **Preconditions:**
    
    - All methods assume that input integers are positive (n > 0) unless noted otherwise.
    - For functions using random sampling (Fermat and Miller–Rabin tests), n should be sufficiently large (typically n > 3)
      to avoid errors in random range generation.
    """

    @staticmethod
    def coprime(a: int, b: int) -> bool:
        """
        Determine whether two integers are coprime (i.e., their greatest common divisor is 1).

        **Preconditions:**
        
        - `a` and `b` are integers (preferably non-negative).
        
        :param a: First integer.
        :param b: Second integer.
        :return: ``True`` if gcd(a, b) == 1, meaning a and b have no common factors other than 1;
                 ``False`` otherwise.
        """
        return math.gcd(a, b) == 1

    @staticmethod
    def extended_euclidean(a: int, b: int) -> tuple[int, int]:
        """
        Compute the coefficients (x, y) for the extended Euclidean algorithm.
        
        This finds integers x and y such that:
        
        .. math:: a \\cdot x + b \\cdot y = \\gcd(a,b)
        
        **Preconditions:**
        
        - `a` and `b` are integers and not both zero.
        
        :param a: First integer.
        :param b: Second integer.
        :return: A tuple (x, y) such that a*x + b*y equals the greatest common divisor of a and b.
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
        """
        Test if a number is prime using a simplified version of the Sieve of Eratosthenes.

        **Mathematical Intuition:**
        
        - The algorithm checks divisibility up to the integer square root of n.
        - It is efficient for moderate-sized inputs but assumes n is a positive integer greater than 1.
        
        **Preconditions:**
        
        - n must be an integer with n > 1.
        
        :param n: An integer greater than 1.
        :return: ``True`` if n is prime; ``False`` otherwise.
        """
        if n % 2 == 0:
            return n == 2
        for v in range(3, math.isqrt(n)+1, 2):
            if n % v == 0:
                return False
        return True

    @staticmethod
    def fermat_primality(n: int) -> bool:
        """
        Test if a number is prime using Fermat's little theorem.

        **Mathematical Intuition:**
        
        - Fermat's little theorem states that for a prime number p and any integer a (1 < a < p),
          :math:`a^{p-1} \\equiv 1 \\pmod{p}`.
        - This method is probabilistic: if the test fails, n is composite; if it passes, n is likely prime
          but may be a Fermat pseudoprime.
          
        **Preconditions:**
        
        - n must be an integer with n > 3.
        - The method assumes n is odd (or equals 2 or 3) because even numbers > 2 are composite.
        - The random base a is chosen in the range [2, math.isqrt(n)]; thus, math.isqrt(n) must be at least 2.
        
        :param n: An integer greater than 3 (or 2 or 3, which are handled explicitly).
        :return: ``True`` if n passes Fermat's test; ``False`` if n fails the test.
        """
        if n == 2 or n == 3:
            return True
        a = random.randint(2, math.isqrt(n))
        if not PrimalityTesting.coprime(a, n):
            return False
        return pow(a, n - 1, n) == 1

    @staticmethod
    def miller_rabin_primality(n: int) -> bool:
        """
        Test if a number is prime using the Miller–Rabin probabilistic primality test.

        **Mathematical Intuition:**
        
        - This algorithm decomposes n-1 into 2^s * d and then uses a randomly chosen base a to test for strong 
          witnesses to compositeness.
        - It is more robust than Fermat's test but still probabilistic.
          
        **Preconditions:**
        
        - n must be an integer with n > 3.
        - The algorithm assumes n is odd (or equals 2 or 3, which are handled explicitly).
        - For n > 3, a random base a is chosen from [2, n-2], so n must be at least 4.
        
        :param n: An odd integer greater than 3 (or 2 or 3, handled explicitly).
        :return: ``True`` if n passes the Miller–Rabin test (likely prime); ``False`` if n fails (composite).
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
    """
    Provides functions related to the prime number theorem (PNT) and approximations of prime counts.

    Several methods are provided for:
    
    - Counting primes exactly using a basic probabilistic test.
    - Estimating the probability that a given number is prime.
    - Summing these probabilities to estimate the total number of primes.
    - Approximating the prime count via integration-based formulas.
    
    **Note:** The functions using Fermat's test may overcount primes due to pseudoprimes.
    """

    @staticmethod
    def π(x: int) -> int:
        """
        Count the number of primes less than or equal to x using a basic Fermat test.
        
        **Mathematical Intuition:**
        
        - This function iterates through all integers from 2 to x,
          testing each for primality with Fermat's test.
        - This is a rudimentary approach and may misidentify some composite numbers as primes.
        
        **Preconditions:**
        
        - x must be an integer greater than or equal to 2.
        
        :param x: An integer, the upper bound for counting primes.
        :return: The count of numbers (from 2 to x) that pass Fermat's primality test.
        """
        count = 0
        for n in range(2, x+1):
            if PrimalityTesting.fermat_primality(n):
                count += 1
        return count

    @staticmethod
    def prob(n: int) -> float:
        """
        Compute the probability that a number n is prime based on a product formula.
        
        **Mathematical Intuition:**
        
        - Uses the heuristic that the probability of a number being prime is
          roughly the product of factors (f-1)/f over primes f up to √n.
        - This is not an exact probability but provides a heuristic estimate.
        
        **Preconditions:**
        
        - n must be an integer greater than or equal to 2.
        
        :param n: An integer to evaluate.
        :return: A floating-point number representing the heuristic probability that n is prime.
        """
        bound = math.floor(n ** 0.5)
        F = [((f-1)/f) for f in range(2, bound+1) if PrimalityTesting.fermat_primality(f)]
        return math.prod(F)

    @staticmethod
    def pnt_prime_count(x: int) -> float:
        """
        Estimate the number of primes up to x using the prime number theorem (PNT) via summation.
        
        **Mathematical Intuition:**
        
        - The PNT implies that the density of primes near x is roughly 1/ln(x).
        - This function sums 1/ln(n) for n from 2 to x as a discrete approximation.
        
        **Preconditions:**
        
        - x must be an integer greater than 1.
        
        :param x: An integer, the upper bound for the approximation.
        :return: A floating-point value representing the approximate number of primes ≤ x.
        """
        acc = 0
        for n in range(2, x+1):
            acc += 1 / math.log(n)
        return acc

    @staticmethod
    def approx_pnt_prime_count(x: int) -> float:
        """
        Approximate the number of primes up to x using the standard PNT formula.
        
        **Mathematical Intuition:**
        
        - The prime number theorem approximates the number of primes ≤ x by x / ln(x).
        - This is a continuous (integral) approximation.
        
        **Preconditions:**
        
        - x must be an integer greater than 1.
        
        :param x: An integer, the upper bound for the approximation.
        :return: The approximate number of primes ≤ x as given by x/ln(x).
        """
        return x / math.log(x)

class Factorization:
    """
    Provides methods for computing the divisors and prime factorization of an integer.
    
    **Preconditions:**
    
    - All methods assume n is a positive integer (n ≥ 1).
    """
    @staticmethod
    def divisors(n: int) -> list[int]:
        """
        Compute all divisors of a given number n.
        
        **Mathematical Intuition:**
        
        - The method computes divisors up to the square root of n, then adds their complementary divisors.
        - The returned list is not guaranteed to be fully sorted.
        
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        
        :param n: A positive integer.
        :return: A list of integers representing all divisors of n.
        """
        res = [a for a in range(1, math.isqrt(n)+1) if n % a ==  0]
        return res + [n // f for f in res[::-1] if f * f != n]

    @staticmethod
    def prime_factors(n: int) -> list[int]:
        """
        Determine the prime factors of a given number n.
        
        **Mathematical Intuition:**
        
        - It computes all divisors (excluding 1) and then filters those that are prime using the
          Sieve of Eratosthenes method.
        - This method may not list repeated factors.
        
        **Preconditions:**
        
        - n must be a positive integer greater than 1.
        
        :param n: A positive integer greater than 1.
        :return: A list of prime numbers that are factors of n.
        """
        return [a for a in Factorization.divisors(n)[1:] if PrimalityTesting.eratosthenes_primality(a)]

    @staticmethod
    def factorize(n: int) -> list[tuple[int, int]]:
        """
        Factorize n into its prime factors and corresponding exponents.
        
        **Mathematical Intuition:**
        
        - For each prime factor, determine the highest exponent such that the prime power divides n.
        - Returns a list of tuples (prime, exponent).
        
        **Preconditions:**
        
        - n must be a positive integer greater than 1.
        
        :param n: A positive integer greater than 1.
        :return: A list of tuples where each tuple consists of a prime factor of n and its exponent.
        """
        def helper(n,a):
            exp = 0
            while n % a == 0:
                n //= a
                exp += 1
            return exp
        return [(f, helper(n, f)) for f in Factorization.prime_factors(n)]

class ArithmeticFunctions:
    """
    Implements several arithmetic functions from number theory, including divisor functions and
    multiplicative functions such as Euler's totient function, Liouville's, and Möbius functions.
    
    **Preconditions:**
    
    - The input n should be a positive integer (n ≥ 1).
    """

    @staticmethod
    def σ(n: int) -> int:
        """
        Compute the sum-of-divisors function σ(n).

        **Mathematical Intuition:**
        
        - σ(n) is defined as the sum of all positive divisors of n.
        
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        
        :param n: A positive integer.
        :return: The sum of all divisors of n.
        """
        return sum(Factorization.divisors(n))
    
    @staticmethod
    def σ_k(n: int, k: int) -> int:
        """
        Compute the generalized divisor function, σₖ(n), which is the sum of the k-th powers of all divisors of n.
        
        **Mathematical Intuition:**
        
        - The function σₖ(n) is defined as:
          
          .. math::
              \\sigma_k(n) = \\sum_{d \\mid n} d^k,
          
          where the sum runs over all positive divisors d of n.
        - This generalizes the standard sum-of-divisors function (σ(n)) which is the case when k = 1.
        
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        - k must be an integer. (k can be zero, positive, or negative; note that for k < 0, divisors are raised to a negative power.)
        
        :param n: A positive integer whose divisors are to be evaluated.
        :param k: An integer exponent applied to each divisor.
        :return: The sum of each divisor of n raised to the power of k.
        """
        return sum(d ** k for d in Factorization.divisors(n))

    @staticmethod
    def Τ(n: int) -> int:
        """
        Compute the divisor-counting function Τ(n).

        **Mathematical Intuition:**
        
        - Τ(n) returns the total number of positive divisors of n.
        
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        
        :param n: A positive integer.
        :return: The count of divisors of n.
        """
        return len(Factorization.divisors(n))

    @staticmethod
    def ω(n: int) -> int:
        """
        Compute the number of distinct prime factors of n.

        **Mathematical Intuition:**
        
        - ω(n) counts each prime factor only once, regardless of multiplicity.
        
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        
        :param n: A positive integer.
        :return: The number of distinct prime factors of n.
        """
        return len(Factorization.prime_factors(n))

    @staticmethod
    def Ω(n: int) -> int:
        """
        Compute the total number of prime factors of n (with multiplicity).

        **Mathematical Intuition:**
        
        - Ω(n) sums the exponents in the prime factorization of n.
        
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        
        :param n: A positive integer.
        :return: The total count of prime factors of n, counting repeated factors.
        """
        return sum([a[1] for a in Factorization.factorize(n)])

    @staticmethod
    def ϕ(n: int):
        """
        Compute Euler's totient function ϕ(n), which counts the number of integers up to n that are coprime with n.

        **Mathematical Intuition:**
        
        - ϕ(n) is calculated by multiplying n by the product of (1 - 1/p) for each distinct prime factor p.
        - This formula is valid for n ≥ 1, with ϕ(1) defined as 1.
        
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        
        :param n: A positive integer.
        :return: The value of Euler's totient function ϕ(n).
        """
        return int(n * math.prod([1-(1/p) for p in Factorization.prime_factors(n)]))

    @staticmethod
    def λ(n: int) -> int:
        """
        Compute the Liouville function λ(n).

        **Mathematical Intuition:**
        
        - λ(n) is defined as (-1)^Ω(n), where Ω(n) is the total number of prime factors of n (with multiplicity).
        - It gives insight into the parity of the total prime factor count.
        
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        
        :param n: A positive integer.
        :return: 1 if Ω(n) is even, -1 if Ω(n) is odd.
        """
        return (-1) ** ArithmeticFunctions.Ω(n)

    @staticmethod
    def μ(n: int) -> int:
        """
        Compute the Möbius function μ(n).

        **Mathematical Intuition:**
        
        - μ(n) is defined as:
          
          - μ(n) = (-1)^ω(n) if n is square-free (i.e., no prime factor is repeated),
          - μ(n) = 0 if n has a squared prime factor.
          
        **Preconditions:**
        
        - n must be a positive integer (n ≥ 1).
        
        :param n: A positive integer.
        :return: The value of the Möbius function: either -1, 0, or 1.
        """
        if ArithmeticFunctions.ω(n) == ArithmeticFunctions.Ω(n):
            return (-1) ** ArithmeticFunctions.Ω(n)
        else:
            return 0

    @staticmethod
    def I(n: int) -> int:
        """
        The indicator function I(n), defined as 1 if n == 1 and 0 otherwise.

        **Mathematical Intuition:**
        
        - This function serves as the identity for Dirichlet convolution in multiplicative number theory.
        
        **Preconditions:**
        
        - n is an integer.
        
        :param n: An integer.
        :return: 1 if n equals 1; 0 otherwise.
        """
        if n == 1:
            return 1
        else:
            return 0

    @staticmethod
    def N(n: int) -> int:
        """
        The identity function N(n) which returns n itself.

        **Mathematical Intuition:**
        
        - In many arithmetic convolutions, N(n) acts as the natural number identity.
        
        **Preconditions:**
        
        - n is an integer.
        
        :param n: An integer.
        :return: The same integer n.
        """
        return n

class ModularAlgebra:
    pass