import sys, tempfile, os, logging, shutil
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from external.adapter import init_hapne, plot_results
from external.hapne.backend.HapNe import IBDLoader, IBDModel


def get_sample_size(vcf_files):
    file = vcf_files[0]
    # Get the number of samples from a vcf.gz file
    # bcftools query -l example/0.vcf.gz | wc -l
    num_samples = subprocess.check_output(
        f"bcftools query -l {file} | wc -l", shell=True
    )
    return int(num_samples)


def read_bins(infiles):
    file = infiles[0]
    return pd.read_csv(file, sep="\t", header=None).values[0:, 0:2].T


def genome_split(region_file):
    return pd.read_csv(region_file, sep="\t")


def load_data(loader, infiles):
    cumulated_data = []
    for filename in infiles:
        data = pd.read_csv(filename, header=None, sep="\t").values[:, 2]
        data = loader.apply_filter_u(data.reshape(-1, 1))
        cumulated_data.append(data.flatten())
    return np.asarray(cumulated_data).T


def compute_poisson_mean(data, region_length):
    p = data.sum(axis=0) / region_length.sum()
    mu = p.reshape((1, -1)) * region_length.reshape(-1, 1)
    return mu


def init_data_handler(infiles, sample_size, region_file, params):
    def lookup(key, default):
        if params.get(key) is not None:
            return params[key]
        return default

    # Skip calling the __init__ method
    loader = IBDLoader.__new__(IBDLoader)
    # Add super class attributes
    loader.genome_split = genome_split(region_file)
    loader.nb_chromosomes = loader.genome_split.shape[0]  # legacy
    loader.u_min = lookup("u_min", loader.get_default_u_min())
    loader.u_max = lookup("u_max", loader.get_default_u_max())
    loader.sample_size = sample_size
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
    # Add IBD specific attributes
    loader.trim = lookup("segment_center_threshold_M", 0.0)
    loader._region_length = loader.genome_split["LENGTH"].values / 100 - 2 * loader.trim
    loader._mu = load_data(loader, infiles).T
    loader.phi2 = np.ones(loader.mu.shape[1])
    loader.adjust_for_count()
    if loader.apply_filter:
        print("Here")
        loader.apply_ibd_region_filter()
        loader.adjust_for_count()
        logging.info(
            f"({loader.genome_split.shape[0] - loader.mu.shape[0]} were discarded)"
        )
    logging.info(
        f"Last bin considered after filtering out small counts: {loader.bins[-1]}"
    )
    logging.info(f"Analyzing {loader.mu.shape[0]} regions ")
    data = loader._mu[loader.selected_regions, :]
    mu = compute_poisson_mean(loader.mu, loader.region_length())
    mu = np.maximum(1.0 / loader.region_length().sum(), mu)
    loader.phi2 = np.var((data - mu) / np.sqrt(mu), axis=0, ddof=1)
    logging.info("Assuming all samples originated from the same generation...")
    loader.time_heterogeneity = None
    return loader


def init_stats_model(loader, sample_size):
    nb_samples = sample_size
    number_haploid_samples = 2 * nb_samples
    nb_pairs = 0.5 * number_haploid_samples * (number_haploid_samples - 1)
    return IBDModel(nb_pairs, loader.region_length().reshape([-1, 1]), loader.phi2)


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    method = "ibd"
    sample_size = get_sample_size(snakemake.input.vcfs)
    infiles = snakemake.input.hists
    region_file = snakemake.input.genome_build
    params = snakemake.params

    data_handler = init_data_handler(infiles, sample_size, region_file, params)
    stats_model = init_stats_model(data_handler, sample_size)
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
        shutil.copy(temp_folder / "hapne.csv", snakemake.output.table)
        hapne.plot_goodness_of_fit(np.median(ne_boot, axis=0), temp_folder)
        shutil.copy(temp_folder / "residuals.png", snakemake.output.residuals)

    with open(str(snakemake.output.summary), "w") as f:
        f.write(hapne.io.summary_message)
    hapne_results = pd.read_csv(snakemake.output.table)
    plot_results(hapne_results, str(snakemake.output.popsize), color="tab:blue")
