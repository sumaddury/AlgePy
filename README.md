# AlgePy: Algebraic Tools in Pure Python

![AlgePy Logo](https://upload.wikimedia.org/wikipedia/commons/thumb/f/fd/Algebra_Proof_Diagram_Inverse.svg/640px-Algebra_Proof_Diagram_Inverse.svg.png)

---

## 1. Introduction

### Module Overview

**AlgePy** is a Python library for discrete algebra and computational number theory. It provides efficient, mathematically rigorous tools for working with number theory and algebraic structures. All algorithms are verified for complexity and correctness.

### Functionality

#### `DiscreteFunctions`
- **`PrimalityTesting`:** Implements deterministic (Sieve of Eratosthenes) and probabilistic (Fermat, Miller–Rabin) methods for testing prime numbers.
- **`PrimeNumberTheorem`:** Provides functions to estimate the density and count of primes.
- **`Factorization`:** Offers algorithms to compute divisors and perform prime factorization.
- **`ArithmeticFunctions`:** Contains number-theoretic functions such as Euler’s totient, Möbius function, Liouville function, and divisor sums.

#### `SingletonStructures`
- **`Z`:**  Represents elements of the ring of integers (ℤ) with overloaded arithmetic and number-theoretic checks.
- **`R`:**  Models real numbers (ℝ) with standard arithmetic operations.
- **`Z_n` & `Z_mod_`:**  Provide modular arithmetic (ℤₙ) functionality, including methods for computing orders, inverses, cyclicity, primitive roots, and the Legendre symbol.

AlgePy also integrates empirical complexity analysis to assess the performance of its algorithms.

---
## 2. Getting Started

### Installation

Install AlgePy directly from GitHub:

```bash
pip install git+https://github.com/sumaddury/AlgePy.git
```

Or install directly from source

```bash
git clone https://github.com/sumaddury/AlgePy.git
cd algepy
pip install .
```

### Quick Usage Example

```python
from AlgePy.DiscreteFunctions import PrimalityTesting, Factorization
from AlgePy.SingletonStructures import Z

# Check if a number is prime using the Sieve of Eratosthenes
n = 29
print(f"{n} is prime: {PrimalityTesting.eratosthenes_primality(n)}")

# Factorize a number
print("Factorization of 360:", Factorization.factorize(360))

# Work with the ring of integers ℤ
a = Z(15)
b = Z(10)
print("a + b =", a + b)
```

---
## 3. Documentation

### Full Documentation: [sumaddury.github.io/AlgePy](https://sumaddury.github.io/AlgePy/)
### Complexity Analysis: [/complexity/complexity_report.pdf](https://github.com/sumaddury/AlgePy/blob/main/complexity/complexity_report.pdf)

---
## 4. Citation

```yaml
cff-version: 1.2.0
title: "AlgePy"
version: "0.0.1"
date-released: "2025-03-31"
authors:
  - given-names: Sucheer
    family-names: Maddury
repository-code: "https://github.com/sumaddury/algepy"
license: "Apache-2.0"
```
