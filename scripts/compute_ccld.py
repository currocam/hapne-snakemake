import logging
from external.hapne.ld import load_and_preprocess_file, _compute_cc_quantities


bed1 = snakemake.input.bed1
bed2 = snakemake.input.bed2
maf = snakemake.params.maf
pseudo_diploid = snakemake.params.pseudo_diploid
nb_points = float(snakemake.params.nb_points)
outfile = str(snakemake.output)
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)

    genotype1, _ = load_and_preprocess_file(bed1, maf)
    genotype2, _ = load_and_preprocess_file(bed2, maf)

    ccld, expected_ccld, bessel_factor, s_corr = _compute_cc_quantities(
        genotype1, genotype2, int(nb_points), pseudo_diploid
    )
    with open(outfile, "w") as f:
        reg1 = snakemake.wildcards.name1
        reg2 = snakemake.wildcards.name2
        f.write(f"{reg1},{reg2},{ccld},{expected_ccld},{bessel_factor},{s_corr}\n")
