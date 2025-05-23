import time
import math
import big_o
import matplotlib.pyplot as plt
import os
import sys

def analyze_complexity(func, 
                         data_gen, 
                         min_n=1000, 
                         max_n=100000, 
                         step=1000, 
                         cluster=1, 
                         n_repeats=50, 
                         plot=False):
    """
    Analyze the time complexity of a given function using the big_o module and additional measurements,
    with support for clustering input sizes to reduce measurement variance.

    The function measures the execution time of the target function (``func``) on inputs generated by
    the provided data generator (``data_gen``) over a range of input sizes. It then groups the measured
    times using a clustering parameter to smooth out noise, and optionally plots the resulting running
    times against the input sizes.

    :param func: Callable
        The target function to analyze. It should accept a single argument, which is generated by
        ``data_gen``.
    :param data_gen: Callable
        A data generator function that accepts an integer n and returns an input of size n.
    :param min_n: int, optional
        The minimum input size for the analysis (default is 1000).
    :param max_n: int, optional
        The maximum input size for the analysis (default is 100000).
    :param step: int, optional
        The step size between successive input sizes (default is 1000).
    :param cluster: int, optional
        The number of successive input sizes to group together (average) to reduce variance.
        If set to 1, no clustering is applied (default is 1).
    :param n_repeats: int, optional
        The number of repetitions to run for each input size when using the big_o estimator (default is 50).
    :param plot: bool, optional
        If True, a plot of the measured execution times (after clustering) is displayed (default is True).

    :returns: tuple
        A tuple containing:
        
        - best_fit: The best fitting complexity model from big_o (e.g., a Linear model).
        - all_fits: A dictionary mapping each tested complexity class to its residual error.
        - ns_clustered: A list of representative input sizes after clustering.
        - times_clustered: A list of the corresponding averaged execution times for these input sizes.
        
    The returned data (``ns_clustered`` and ``times_clustered``) can be used to generate custom plots and
    reports. The output of the big_o estimator (``best_fit`` and ``all_fits``) provides detailed information
    on the estimated asymptotic behavior of the function.
    """
    best_fit, all_fits = big_o.big_o(func, data_gen, n_repeats=n_repeats, min_n=min_n, max_n=max_n)
    
    ns = list(range(min_n, max_n + 1, step))
    raw_times = []
    
    for n in ns:
        start_time = time.time()
        func(data_gen(n))
        end_time = time.time()
        raw_times.append(end_time - start_time)
    
    if cluster > 1:
        ns_clustered = []
        times_clustered = []
        for i in range(0, len(ns), cluster):
            cluster_ns = ns[i:i+cluster]
            cluster_times = raw_times[i:i+cluster]
            representative_n = sum(cluster_ns) / len(cluster_ns)
            avg_time = sum(cluster_times) / len(cluster_times)
            ns_clustered.append(representative_n)
            times_clustered.append(avg_time)
    else:
        ns_clustered = ns
        times_clustered = raw_times

    if plot:
        plt.figure(figsize=(8, 5))
        plt.plot(ns_clustered, times_clustered, 'bo-', label='Measured Time')
        plt.xlabel("Input Size (n)")
        plt.ylabel("Execution Time (sec)")
        plt.title(f"Execution Time vs. Input Size for {func.__name__}")
        plt.legend()
        plt.grid(True)
        plt.show()
    
    return best_fit, all_fits, ns_clustered, times_clustered

def safe_data_gen(data_gen, force_int=False):
    def wrapper(n):
        value = data_gen(n)
        if force_int:
            try:
                return int(value)
            except Exception:
                return value
        return value
    return wrapper

def save_plot(ns, times, title, filename):
    plt.figure(figsize=(8, 5))
    plt.plot(ns, times, 'bo-', label='Measured Time')
    plt.xlabel("Input Size (n)")
    plt.ylabel("Execution Time (sec)")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(filename)
    plt.close()

def generate_report(target_name, best_fit, all_fits, plot_filename):
    md = []
    md.append(f"## Complexity Analysis for `{target_name}`\n\n")
    md.append("**Best Fit Complexity:**\n")
    md.append(f"`{best_fit}`\n\n")
    md.append("**Detailed Fit Residuals:**\n")
    md.append("| Complexity Class      | Residual |\n")
    md.append("|-----------------------|----------|\n")
    for comp, res in all_fits.items():
        md.append(f"| {str(comp):<21} | {res:.2G} |\n")
    md.append("\n**Execution Time vs. Input Size Plot:**\n")
    md.append(r"\begin{center}" + "\n")
    md.append(f"\\includegraphics[width=0.8\\textwidth]{{{plot_filename}}}\n")
    md.append(r"\end{center}" + "\n")
    md.append("\n---\n")
    return "".join(md)