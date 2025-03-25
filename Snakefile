import pandas as pd
from itertools import combinations


# Configuration of the pipeline
configfile: "config.yaml"


data_dir = config["data_dir"]
# Use glob_wildcards to find all VCF files in the directory
# This will match files like '/path/to/data/file1.vcf.gz', '/path/to/data/file2.vcf.gz', etc.
(names,) = glob_wildcards(f"{data_dir}{{name}}.vcf.gz")
# Compute unique pairs of regions from the names vector
unique_pairs = [f"{x1}_{x2}" for x1, x2 in combinations(names, 2)]
nb_regions = len(names)
assert nb_regions * (nb_regions - 1) // 2 == len(unique_pairs)
if config["method"] == "ld":
    all_files = ["steps/hapne.ccld", expand("steps/{name}.r2", name=names)]
else:
    raise NotImplementedError("Only the 'ld' method is supported")


rule all:
    input:
        all_files,


# LD-based pipeline
rule split_vcf:
    input:
        vcf=f"{data_dir}{{name}}.vcf.gz",
        map=f"{data_dir}{{name}}.map",
    output:
        "steps/{name}.bed",
        "steps/{name}.bim",
        "steps/{name}.fam",
    log:
        "logs/plink_{name}.log",
    threads: 1
    shadow:
        "minimal"
    shell:
        """
        plink --vcf {input.vcf} \
            --make-bed \
            --chr {wildcards.name} \
            --cm-map {input.map} {wildcards.name} \
            --out steps/{wildcards.name} \
            --threads {threads} \
            --memory 2048 \
            --maf 0.249 \
            --snps-only \
            --geno 0.8 \
            > {log}
        """


rule compute_ld:
    input:
        bed="steps/{name}.bed",
    output:
        "steps/{name}.r2",
    params:
        maf=0.25,
        pseudo_diploid=False,
        bin_file=config.get("bin_file"),
    log:
        "logs/ld_{name}.log",
    threads: 1
    script:
        "scripts/compute_ld.py"


rule compute_ccld_pair:
    input:
        bed1="steps/{name1}.bed",
        bed2="steps/{name2}.bed",
    output:
        temp("steps/{name1}_{name2}.r2"),
    log:
        "logs/ld_{name1}_{name2}.log",
    params:
        maf=config.get("maf", 0.25),
        pseudo_diploid=config.get("pseudo_diploid", False),
        nb_points=config.get("nb_points", int(1e6)),
    threads: 1
    script:
        "scripts/compute_ccld.py"


rule compute_ccld:
    input:
        expand("steps/{name}.r2", name=unique_pairs),
    output:
        "steps/hapne.ccld",
    shell:
        """
        echo REGION1, REGION2, CCLD, CCLD_H0, BESSEL_FACTOR, S_CORR > {output}
        cat {input} >> {output}
        """
