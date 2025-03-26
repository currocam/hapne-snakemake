import pandas as pd
from itertools import combinations


# Configuration of the pipeline
configfile: "config.yaml"


data_dir = config["data_dir"]
out_dir = config["out_dir"]
# Use glob_wildcards to find all VCF files in the directory
# This will match files like '/path/to/data/file1.vcf.gz', '/path/to/data/file2.vcf.gz', etc.
(names,) = glob_wildcards(f"{data_dir}{{name}}.vcf.gz")
# Compute unique pairs of regions from the names vector
unique_pairs = [f"{x1}_{x2}" for x1, x2 in combinations(names, 2)]
nb_regions = len(names)
assert nb_regions * (nb_regions - 1) // 2 == len(unique_pairs)
if config["method"] == "ld":
    all_files = [
        f"{out_dir}ld_hapne_estimate.csv",
        f"{out_dir}ld_hapne_summary.txt",
        f"{out_dir}ld_hapne_residuals.png",
        f"{out_dir}ld_hapne_pop_trajectory.png",
    ]
else:
    raise NotImplementedError("Only the 'ld' method is supported")


rule all:
    input:
        all_files,


# LD-based pipeline
rule convert_plink:
    input:
        vcf=f"{data_dir}{{name}}.vcf.gz",
        map=f"{data_dir}{{name}}.map",
    output:
        "steps/{name}.bed",
        "steps/{name}.bim",
        "steps/{name}.fam",
    log:
        "logs/plink/{name}.log",
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
        "logs/ld/{name}.log",
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
        "logs/ld/{name1}_{name2}.log",
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


rule fit_hapne_ld:
    input:
        infiles=expand("steps/{name}.r2", name=names),
        fam_files=expand("steps/{name}.fam", name=names),
        bias="steps/hapne.ccld",
    params:
        apply_filter=config.get("apply_filter", True),
        pseudo_diploid=config.get("pseudo_diploid", False),
        nb_individuals=config.get("nb_individuals"),
        u_min=config.get("u_min"),
        u_max=config.get("u_max"),
        filter_tol=config.get("filter_tol"),
        sigma2=config.get("sigma2"),
        u_quantile=config.get("u_quantile"),
        dt_min=config.get("dt_min"),
        dt_max=config.get("dt_max"),
        t_max=config.get("t_max"),
        nb_parameters=config.get("nb_parameters"),
        model=config.get("mode"),
    log:
        "logs/ld_hapne.log",
    output:
        table=f"{out_dir}ld_hapne_estimate.csv",
        summary=f"{out_dir}ld_hapne_summary.txt",
        residuals=f"{out_dir}ld_hapne_residuals.png",
        popsize=f"{out_dir}ld_hapne_pop_trajectory.png",
    script:
        "scripts/fit_hapne_ld.py"
