from external.hapne.backend.HapNe import HapNe, LDLoader, LDModel
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import tempfile
from configparser import ConfigParser
import os
import numpy as np
import logging
from pathlib import Path
from scipy.stats import ttest_1samp


# The neccesary methods have been re-implemented here
def load_data(infiles):
    cumulated_data = [pd.read_csv(file)["R2"].values for file in infiles]
    return np.asarray(cumulated_data).T


def get_sample_size(fam_files):
    pseudo_diploid = lookup("pseudo_diploid", False)
    if not pseudo_diploid:
        logging.warning(
            "Assuming diploid. If you are analysing aDNA, set the flag to True."
        )
    try:
        famfile = fam_files[0]
        fam_file = pd.read_csv(famfile, header=None, sep="\t")
        nb_individuals = fam_file.shape[0]
    except FileNotFoundError:
        nb_individuals = lookup("nb_individuals", 0)
        if nb_individuals == 0:
            raise ValueError("Number of individuals not provided")
        nb_individuals = nb_individuals
    return (2 - pseudo_diploid) * nb_individuals


def read_bins(infiles):
    file = infiles[0]
    bins = pd.read_csv(file).values[:, 1:3].T
    return bins


def genome_split(infiles):
    # Basename of the files
    names = [os.path.basename(x) for x in infiles]
    return pd.DataFrame(
        {
            "CHR": names,
            "FROM_BP": np.nan,
            "TO_BP": np.nan,
            "NAME": names,
            "LENGTH": np.nan,
        }
    )


def init_loader(infiles, fam_files, bias_filename):
    # Skip calling the __init__ method
    loader = LDLoader.__new__(LDLoader)
    loader.genome_split = genome_split(infiles)
    loader.nb_chromosomes = len(infiles)
    loader.u_min = lookup("u_min", loader.get_default_u_min())
    loader.u_max = lookup("u_max", loader.get_default_u_max())
    loader.sample_size = get_sample_size(fam_files)
    message = f"Analysing {loader.sample_size} "
    message += "haplotypes"
    logging.info(message)
    loader._bins = read_bins(infiles)
    bin_from = loader.apply_filter_u(loader._bins[0, :].reshape(-1, 1))
    bin_to = loader.apply_filter_u(loader._bins[1, :].reshape(-1, 1))
    loader.bins = np.append(bin_from.reshape(-1), bin_to[-1, 0]).reshape(-1)
    loader._mu = None
    loader.selected_regions = np.arange(loader.genome_split.shape[0])
    loader.apply_filter = lookup("apply_filter", True)
    loader.summary_message = ""
    # Read the data (from .load_data method)
    loader.ld = load_data(infiles)
    # Logic about the bias file
    data = pd.read_csv(bias_filename)
    # Original code rely on the specific format of the regions chr1, chr2, etc.
    # Here, I assume each region is from a different chromosome
    logging.warning("Assuming each region is from a different chromosome")
    # Select regions where the first 4 character do not match
    # region_1 = data.values[:, 0]
    # region_2 = data.values[:, 1]
    # valid_rows = np.array([region_1[i][:5] != region_2[i][:5] for i in range(len(region_1))])
    valid_rows = np.ones(data.shape[0], dtype=bool)
    loader.ccld = data.values[valid_rows, 2].mean()
    loader.bessel_cor = data.values[valid_rows, 4].mean()
    _, loader.pval_no_admixture = ttest_1samp(
        data.values[valid_rows, 2].astype(float)
        - data.values[valid_rows, 3].astype(float),
        0,
    )
    loader.summary_message += f"CCLD: {loader.ccld:.5f}. \n"
    loader.summary_message += f"The p-value associated with H0 = no structure is {loader.pval_no_admixture:.3f}.\n"
    loader.summary_message += (
        "If H0 is rejected, contractions in the recent past "
        "might reflect structure instead of reduced population size."
    )
    logging.warning(f"CCLD: {loader.ccld:.5f}.")
    logging.warning(
        f"The p-value associated with H0 = no structure is {loader.pval_no_admixture:.3f}."
    )
    logging.warning(
        (
            "If H0 is rejected, contractions in the recent past might reflect "
            "structure instead of reduced population size."
        )
    )
    loader._mu = loader.ld_to_p_ibd(loader.apply_filter_u(loader.ld)).T
    # Combine bins if we are dealing with a small genotype or high missingness
    # loader.merge_bins_if_required()
    # Discard suspicious Regions
    loader.filter_tolerance = lookup("filter_tol", 5e-4)
    if loader.apply_filter:
        loader.apply_region_filter()
    logging.info(f"Analyzing {loader.mu.shape[0]} regions ")

    loader.phi2 = np.var(loader.mu, axis=0).reshape(-1)
    loader.sigma = np.cov(loader.mu.T)

    logging.info("Assuming all samples originated from the same generation...")
    loader.time_heterogeneity = None
    return loader


def init_hapne(infiles, fam_files, bias_file, params):
    # Skip calling the __init__ method
    hapne = HapNe.__new__(HapNe)
    # The object expects a very specific ConfigParser object
    hapne.config = ConfigParser()
    hapne.config.add_section("CONFIG")
    # Iterate over params and add them to the config
    for key, value in params.items():
        hapne.config.set("CONFIG", key, str(value))
    hapne.method = "ld"
    hapne.io = init_loader(infiles, fam_files, bias_file)
    hapne.stats_model = LDModel(sigma=hapne.io.sigma)
    # Parameter related to the model
    nb_points = 7
    start = -5
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
    hapne.signal_threshold = 6.635 if hapne.method == "ld" else 9.210
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


# Get data from snakemake object
infiles = snakemake.input.infiles
fam_files = snakemake.input.fam_files
bias_filename = snakemake.input.bias
params = snakemake.params


# Global variable with configuration
def lookup(key, default):
    if params.get(key) is not None:
        return params[key]
    return default


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    hapne = init_hapne(infiles, fam_files, bias_filename, params)
    if hapne.mode == "regularised":
        ne_boot = hapne.fit_regularized()
    elif hapne.mode == "fixed":
        ne_boot = hapne.fit_fixed()
    else:
        raise ValueError("Unknown mode")

    # Create temp folder and save results
    with tempfile.TemporaryDirectory() as temp_folder:
        temp_folder = Path(temp_folder)
        hapne.config.set("CONFIG", "output_folder", str(temp_folder))
        hapne.save_hapne(ne_boot, temp_folder)
        os.rename(temp_folder / "hapne.csv", snakemake.output.table)
        hapne.plot_goodness_of_fit(np.median(ne_boot, axis=0), temp_folder)
        os.rename(temp_folder / "residuals.png", snakemake.output.residuals)
    with open(str(snakemake.output.summary), "w") as f:
        f.write(hapne.io.summary_message)
    hapne_results = pd.read_csv(snakemake.output.table)
    plot_results(hapne_results, str(snakemake.output.popsize), color="tab:blue")
