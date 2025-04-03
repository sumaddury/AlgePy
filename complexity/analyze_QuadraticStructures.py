import os
import time
import matplotlib.pyplot as plt
import sys
import subprocess
from Analyze import analyze_complexity, safe_data_gen, save_plot, generate_report

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'algepy')))

REPORTS_DIR = os.path.join(os.path.dirname(__file__), "reports/complexityQuadraticStructures")
os.makedirs(REPORTS_DIR, exist_ok=True)

# Define analysis targets.
# Each target is a dictionary with keys:
#   - 'name': a descriptive name,
#   - 'func': the function to analyze,
#   - 'data_gen': a data generator that accepts an integer n and returns an input,
#   - 'min_n', 'max_n', 'step', 'cluster', 'n_repeats': analysis parameters,
#   - 'force_int': if True, the generated value is cast to int.
from SingletonStructures import Z, R
from QuadraticStructures import C, Q, QuadInt, QuadIntRing, QuadRat, QuadRatField
ring = QuadIntRing(2, force_ufd=True)
ring_imag = QuadIntRing(-2, force_ufd=True)
field = QuadRatField(2)


targets = [
    {
        "name": "quadint_gcd",
        "func": lambda pair: ring.gcd(*pair),
        "data_gen": lambda n: (ring(int(n), 0), ring(int(n) + 14, 0)),
        "min_n": 10,
        "max_n": 10000,
        "step": 1,
        "cluster": 50,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "quadint_exponentiation",
        "func": lambda z: z ** Z(5),
        "data_gen": lambda n: ring(int(n), 0),
        "min_n": 1,
        "max_n": 10000,
        "step": 1,
        "cluster": 50,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "quadint_floor_division",
        "func": lambda pair: ring.internal_div(*pair),
        "data_gen": lambda n: (ring(int(n), 0), ring(int(n) // 2 if int(n) // 2 > 0 else 1, 0)),
        "min_n": 10,
        "max_n": 10000,
        "step": 1,
        "cluster": 50,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "quadint_embedding",
        "func": lambda z: ring.embed(z),
        "data_gen": lambda n: ring(int(n), int(n) + 1),
        "min_n": 10,
        "max_n": 10000,
        "step": 1,
        "cluster": 50,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "factorization_standard",
        "func": lambda z: ring.factorize(z),
        "data_gen": lambda n: ring(int(n), 0),  
        "min_n": 2,
        "max_n": 10000,
        "step": 5,
        "cluster": 20,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "factorization_imaginary",
        "func": lambda z: ring_imag.factorize(z),
        "data_gen": lambda n: ring_imag(int(n), 0),  
        "min_n": 2,
        "max_n": 10000,
        "step": 5,
        "cluster": 20,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "normalize_standard",
        "func": lambda z: z.normalize(),
        "data_gen": lambda n: ring(int(n), int(n) // 3 if int(n) // 3 != 0 else 1),
        "min_n": 10,
        "max_n": 10000,
        "step": 5,
        "cluster": 20,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "normalize_imaginary",
        "func": lambda z: z.normalize(),
        "data_gen": lambda n: ring_imag(int(n), int(n) // 3 if int(n) // 3 != 0 else 1),
        "min_n": 10,
        "max_n": 10000,
        "step": 5,
        "cluster": 20,
        "n_repeats": 10,
        "force_int": True
    },

    # ------------------ Quadratic Rational Field (d=2) ------------------
    {
        "name": "quadrat_exponentiation",
        "func": lambda q: q ** Z(5),
        "data_gen": lambda n: field(Q(int(n), 2), Q(int(n) + 2, 2)),
        "min_n": 1,
        "max_n": 10000,
        "step": 1,
        "cluster": 50,
        "n_repeats": 10,
        "force_int": True
    },
    {
        "name": "quadrat_embedding",
        "func": lambda q: field.embed(q),
        "data_gen": lambda n: field(Q(int(n), 2), Q(int(n) + 2, 2)),
        "min_n": 1,
        "max_n": 10000,
        "step": 1,
        "cluster": 50,
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
  - \\graphicspath{{reports/complexityQuadraticStructures/}}
---
"""

aggregated_report = header + "\n" + "# Quadratic Structures Complexity Analysis Report\n\n" + "\n".join(all_reports)
report_filepath = os.path.join(REPORTS_DIR, "quadratic_structures_complexity_report.md")
with open(report_filepath, "w") as f:
    f.write(aggregated_report)

print(f"Aggregated complexity report saved to: {report_filepath}")

pdf_output = os.path.join(os.path.dirname(__file__), "quadratic_structures_complexity_report.pdf")

subprocess.run([
    "pandoc", report_filepath, "-o", pdf_output,
    "--pdf-engine=tectonic"
], check=True)

print(f"PDF report saved to: {pdf_output}")