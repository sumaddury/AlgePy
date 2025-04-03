import os
import time
import matplotlib.pyplot as plt
import sys
import subprocess
from Analyze import analyze_complexity, safe_data_gen, save_plot, generate_report

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'algepy')))

REPORTS_DIR = os.path.join(os.path.dirname(__file__), "reports/complexitySingletonStructures")
os.makedirs(REPORTS_DIR, exist_ok=True)

# Define analysis targets.
# Each target is a dictionary with keys:
#   - 'name': a descriptive name,
#   - 'func': the function to analyze,
#   - 'data_gen': a data generator that accepts an integer n and returns an input,
#   - 'min_n', 'max_n', 'step', 'cluster', 'n_repeats': analysis parameters,
#   - 'force_int': if True, the generated value is cast to int.
from .SingletonStructures import Z, R, Z_n, Z_mod_

targets = [
    {
        "name": "Z_n_order",
        "func": lambda z: z.order(),
        "data_gen": lambda n: Z_n(int(n), 101),
        "min_n": 1,
        "max_n": 100,
        "step": 1,
        "cluster": 1,
        "n_repeats": 10,
        "force_int": False
    },
    {
        "name": "Z_n_inverse",
        "func": lambda z: z.inverse(),
        "data_gen": lambda n: Z_n(int(n), 101),
        "min_n": 1,
        "max_n": 100,
        "step": 1,
        "cluster": 1,
        "n_repeats": 10,
        "force_int": False
    },
    {
        "name": "Z_n_legendre",
        "func": lambda x: Z_mod_(125371).legendre(Z(int(x))),
        "data_gen": lambda n: n,
        "min_n": 1,
        "max_n": 125370,
        "step": 1,
        "cluster": 1,
        "n_repeats": 10,
        "force_int": True
    },
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
  - \\graphicspath{{reports/complexitySingletonStructures/}}
---
"""

aggregated_report = header + "\n" + "# Singleton Structures Complexity Analysis Report\n\n" + "\n".join(all_reports)
report_filepath = os.path.join(REPORTS_DIR, "singleton_structures_complexity_report.md")
with open(report_filepath, "w") as f:
    f.write(aggregated_report)

print(f"Aggregated complexity report saved to: {report_filepath}")

pdf_output = os.path.join(os.path.dirname(__file__), "singleton_structures_complexity_report.pdf")

subprocess.run([
    "pandoc", report_filepath, "-o", pdf_output,
    "--pdf-engine=tectonic"
], check=True)

print(f"PDF report saved to: {pdf_output}")