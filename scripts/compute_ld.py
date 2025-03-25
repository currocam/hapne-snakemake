import numpy as np
import pandas_plink as pdp
import logging
import pandas as pd
from external.hapne.utils import get_bins
from external.hapne.ld import (
    load_and_preprocess_file,
    _compute_r2,
)

bed_file = snakemake.input.bed
maf = snakemake.params.maf
pseudo_diploid = snakemake.params.pseudo_diploid
outfile = str(snakemake.output)

logger = logging.getLogger(__name__)
logging.basicConfig(filename=str(snakemake.log), encoding="utf-8", level=logging.INFO)

genotype, chr_map = load_and_preprocess_file(bed_file, maf)
average_missing_prop = np.mean(np.isnan(genotype))
sample_size = genotype.shape[0] * (2 - pseudo_diploid)
if sample_size * (1 - average_missing_prop) ** 2 < 6:
    raise ValueError(
        f"File {bed_file} has too many missing values or too few individuals."
    )
else:
    logging.info(
        f"Computing LD for file {bed_file} with {genotype.shape[0]} individuals and {genotype.shape[1]} SNPs."
    )
bins = snakemake.params.get("bin_file")
if bins is None:
    bins = get_bins()
else:
    # Read from a file instead
    bins = np.loadtxt(bins)

ld, weights = _compute_r2(genotype, chr_map, bins)

to_be_saved = pd.DataFrame(columns=["BIN_FROM[M]", "BIN_TO[M]", "R2", "WEIGHT"])
to_be_saved["BIN_FROM[M]"] = bins[:, 0]
to_be_saved["BIN_TO[M]"] = bins[:, 1]
to_be_saved["R2"] = ld
to_be_saved["WEIGHT"] = weights
to_be_saved.to_csv(outfile)
