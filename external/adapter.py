from external.hapne.backend.HapNe import HapNe
import numpy as np
import scipy
from configparser import ConfigParser
import logging
import matplotlib.pyplot as plt

def init_hapne(method, params, data_handler, stats_model):
    def lookup(key, default):
        if params.get(key) is not None:
            return params[key]
        return default
    # Skip calling the __init__ method
    hapne = HapNe.__new__(HapNe)
    # The object expects a very specific ConfigParser object
    hapne.config = ConfigParser()
    hapne.config.add_section("CONFIG")
    # Iterate over params and add them to the config
    for key, value in params.items():
        if value is not None:
            hapne.config.set("CONFIG", key, str(value))
    hapne.method = method
    # Method specific instances
    hapne.io = data_handler
    hapne.stats_model = stats_model
    # Parameter related to the model
    if hapne.method == "ld":
        nb_points = 7
        start = -5
    elif hapne.method == "ibd":
        nb_points = 8
        start = -6
    else:
        raise ValueError("Unknown method")
    if lookup("sigma2", None) is None:
        hapne.parameter_grid = np.array(
            [10.0 ** (start + ii) for ii in range(0, nb_points)]
        )
    else:
        logging.info("Using user-provided sigma2")
        sigma2 = lookup("sigma2", None)
        hapne.parameter_grid = np.array([sigma2])
    # This seems to ignore the user default, but it was done in the original code
    hapne.u_min = hapne.get_default_u_min()
    hapne.u_quantile = lookup("u_quantile", hapne.u_min)
    hapne.dt_min = max(0.25, lookup("dt_min", 1))
    hapne.dt_max = lookup("dt_max", 5000)
    hapne.t_max = lookup("t_max", 125)
    hapne.buffer_time = hapne.t_max - 25  # the times averages the deep time effect
    # Parameters related to the fitting procedure
    hapne.mode = lookup("mode", "regularised")
    pseudo_diploid = lookup("pseudo_diploid", False)
    if hapne.mode == "regularised":
        default_params = 16 if pseudo_diploid else 21
    else:
        default_params = 50
    hapne.nb_parameters = lookup("nb_parameters", default_params)
    hapne.random_restarts = 3
    hapne.eps = 1e-4
    hapne.ftol = 1e-3
    hapne.n_min = 1e2
    hapne.n_max = 1e9
    if hapne.method == "ld":
        hapne.signal_threshold = 6.635
    elif hapne.method == "ibd":
        hapne.signal_threshold = 9.210
    else:
        raise ValueError("Unknown method")
    hapne.nb_bootstraps = lookup("nb_bootstraps", 100)
    # Initialisations
    hapne.times = scipy.stats.erlang.ppf(
        np.linspace(0, 1, hapne.nb_parameters + 1)[:-1],
        a=2,
        scale=1.0 / (2 * hapne.u_quantile),
    )
    hapne.init_times()
    return hapne

def plot_results(hapne_results, outfile, color="tab:blue"):
    time = hapne_results["TIME"]
    xlabel = "Time (gen.)"
    fig, ax = plt.subplots(figsize=(5, 2.3))
    ax.plot(time, hapne_results["Q0.5"], color=color, label="")
    ax.fill_between(
        time, hapne_results["Q0.25"], hapne_results["Q0.75"], color=color, alpha=0.5
    )
    ax.fill_between(
        time, hapne_results["Q0.025"], hapne_results["Q0.975"], color=color, alpha=0.25
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("$N_e$")

    ax.set_yscale("log")
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    plt.savefig(outfile, bbox_inches="tight")
    plt.close()
