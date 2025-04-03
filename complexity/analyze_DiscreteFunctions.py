import os
import time
import matplotlib.pyplot as plt
import sys
import subprocess
from Analyze import analyze_complexity, safe_data_gen, save_plot, generate_report

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

REPORTS_DIR = os.path.join(os.path.dirname(__file__), "reports/complexityDiscreteFunctions")
os.makedirs(REPORTS_DIR, exist_ok=True)

# Define analysis targets.
# Each target is a dictionary with keys:
#   - 'name': a descriptive name,
#   - 'func': the function to analyze,
#   - 'data_gen': a data generator that accepts an integer n and returns an input,
#   - 'min_n', 'max_n', 'step', 'cluster', 'n_repeats': analysis parameters,
#   - 'force_int': if True, the generated value is cast to int.
from DiscreteFunctions import PrimalityTesting, PrimeNumberTheorem, Factorization, ArithmeticFunctions

targets = [
    # --- PrimalityTesting methods ---
    {
        "name": "extended_euclidean",
        "func": lambda x: PrimalityTesting.extended_euclidean(*x),
        "data_gen": lambda n: (int(n), n // 2),
        "min_n": 100,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "eratosthenes_primality",
        "func": PrimalityTesting.eratosthenes_primality,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 20,
        "force_int": True
    },
    {
        "name": "fermat_primality",
        "func": PrimalityTesting.fermat_primality,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 50,
        "force_int": True
    },
    {
        "name": "miller_rabin_primality",
        "func": PrimalityTesting.miller_rabin_primality,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 50,
        "force_int": True
    },
    # --- PrimeNumberTheorem methods ---
    {
        "name": "pi_function",
        "func": PrimeNumberTheorem.π,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 5000,
        "cluster": 1,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "prob_function",
        "func": PrimeNumberTheorem.prob,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 5000,
        "cluster": 1,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "pnt_prime_count",
        "func": PrimeNumberTheorem.pnt_prime_count,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 5000,
        "cluster": 1,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "approx_pnt_prime_count",
        "func": PrimeNumberTheorem.approx_pnt_prime_count,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 5000,
        "cluster": 1,
        "n_repeats": 10,
        "force_int": True
    },
    # --- Factorization methods ---
    {
        "name": "divisors",
        "func": Factorization.divisors,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "prime_factors",
        "func": Factorization.prime_factors,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "factorize",
        "func": Factorization.factorize,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    # --- ArithmeticFunctions methods ---
    {
        "name": "sigma",
        "func": ArithmeticFunctions.σ,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "sigma_k",
        "func": lambda n: ArithmeticFunctions.σ_k(n, 2),
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "T_divisor_count",
        "func": ArithmeticFunctions.Τ,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "omega_distinct",
        "func": ArithmeticFunctions.ω,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "Omega_total",
        "func": ArithmeticFunctions.Ω,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "phi_totient",
        "func": ArithmeticFunctions.ϕ,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "lambda_liouville",
        "func": ArithmeticFunctions.λ,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "mu_mobius",
        "func": ArithmeticFunctions.μ,
        "data_gen": lambda n: n,
        "min_n": 1000,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "I_indicator",
        "func": ArithmeticFunctions.I,
        "data_gen": lambda n: n,
        "min_n": 10,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "N_identity",
        "func": ArithmeticFunctions.N,
        "data_gen": lambda n: n,
        "min_n": 10,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "is_perfect",
        "func": ArithmeticFunctions.is_perfect,
        "data_gen": lambda n: n,
        "min_n": 10,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "is_square_free",
        "func": ArithmeticFunctions.is_square_free,
        "data_gen": lambda n: n,
        "min_n": 10,
        "max_n": 200000,
        "step": 10,
        "cluster": 100,
        "n_repeats": 10,
        "force_int": True
    }
]

all_reports = []

for target in targets:
    name = target["name"]
    func = target["func"]
    data_gen_wrapped = safe_data_gen(target["data_gen"], force_int=target.get("force_int", False))
    min_n = target["min_n"]
    max_n = target["max_n"]
    step = target["step"]
    cluster = target["cluster"]
    n_repeats = target["n_repeats"]

    print(name)

    best_fit, all_fits, ns_clustered, times_clustered = analyze_complexity(
        func, data_gen_wrapped,
        min_n=min_n, max_n=max_n,
        step=step, cluster=cluster,
        n_repeats=n_repeats,
        plot=False
    )   
    
    
    plot_filename = os.path.join(REPORTS_DIR, f"{name}_plot.png")
    plot_title = f"Execution Time vs. Input Size for {name}"
    save_plot(ns_clustered, times_clustered, plot_title, plot_filename)
    
    report_section = generate_report(name, best_fit, all_fits, os.path.basename(plot_filename))
    all_reports.append(report_section)

header = """---
header-includes:
  - \\usepackage{graphicx}
  - \\usepackage{float}
  - \\graphicspath{{reports/complexityDiscreteFunctions/}}
---
"""

aggregated_report = header + "\n" + "# Discrete Functions Complexity Analysis Report\n\n" + "\n".join(all_reports)
report_filepath = os.path.join(REPORTS_DIR, "discrete_functions_complexity_report.md")
with open(report_filepath, "w") as f:
    f.write(aggregated_report)

print(f"Aggregated complexity report saved to: {report_filepath}")

pdf_output = os.path.join(os.path.dirname(__file__), "discrete_functions_complexity_report.pdf")

subprocess.run([
    "pandoc", report_filepath, "-o", pdf_output,
    "--pdf-engine=tectonic"
], check=True)

print(f"PDF report saved to: {pdf_output}")