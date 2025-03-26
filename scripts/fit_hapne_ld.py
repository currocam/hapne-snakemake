from external.hapne.backend.HapNe import HapNe, LDLoader, LDModel
from external.adapter import init_hapne, plot_results
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


def get_sample_size(params, fam_files):
    def lookup(key, default):
        if params.get(key) is not None:
            return params[key]
        return default

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


def init_data_handler(infiles, fam_files, bias_filename, params):
    def lookup(key, default):
        if params.get(key) is not None:
            return params[key]
        return default

    # Skip calling the __init__ method
    loader = LDLoader.__new__(LDLoader)
    loader.genome_split = genome_split(infiles)
    loader.nb_chromosomes = len(infiles)
    loader.u_min = lookup("u_min", loader.get_default_u_min())
    loader.u_max = lookup("u_max", loader.get_default_u_max())
    loader.sample_size = get_sample_size(params, fam_files)
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


def init_stats_model(loader):
    return LDModel(sigma=loader.sigma)


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    method = "ld"
    infiles = snakemake.input.infiles
    fam_files = snakemake.input.fam_files
    bias_filename = snakemake.input.bias
    params = snakemake.params

    data_handler = init_data_handler(infiles, fam_files, bias_filename, params)
    stats_model = init_stats_model(data_handler)
    hapne = init_hapne(method, params, data_handler, stats_model)
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
