import pandas as pd
from itertools import combinations

# Configuration of the pipeline
# configfile: "config.yaml"
data_dir = config["data_dir"]
out_dir = config["out_dir"]
# Use glob_wildcards to find all VCF files in the directory
# This will match files like '/path/to/data/file1.vcf.gz', '/path/to/data/file2.vcf.gz', etc.
(names,) = glob_wildcards(f"{data_dir}{{name}}.vcf.gz")
# Compute unique pairs of regions from the names vector
unique_pairs = [f"{x1}_{x2}" for x1, x2 in combinations(names, 2)]
nb_regions = len(names)
assert nb_regions * (nb_regions - 1) // 2 == len(unique_pairs)
all_files = [
    f"{out_dir}{config["method"]}_hapne_estimate.csv",
    f"{out_dir}{config["method"]}_hapne_summary.txt",
    f"{out_dir}{config["method"]}_hapne_residuals.png",
    f"{out_dir}{config["method"]}_hapne_pop_trajectory.png",
]
suffix_map = config.get("map_file_suffix", ".map")


rule all:
    input:
        all_files,


# LD-based pipeline
rule convert_plink:
    input:
        vcf=f"{data_dir}{{name}}.vcf.gz",
        map=f"{data_dir}{{name}}{suffix_map}",
    output:
        "steps/{name}.bed",
        "steps/{name}.bim",
        "steps/{name}.fam",
    log:
        "logs/plink/{name}.log",
    threads: 1
    conda:
        "conda_environment.yaml"
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
    conda:
        "conda_environment.yaml"
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
    conda:
        "conda_environment.yaml"
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
        nb_bootstraps=config.get("nb_bootstraps"),
    log:
        "logs/ld_hapne.log",
    conda:
        "conda_environment.yaml"
    output:
        table=f"{out_dir}ld_hapne_estimate.csv",
        summary=f"{out_dir}ld_hapne_summary.txt",
        residuals=f"{out_dir}ld_hapne_residuals.png",
        popsize=f"{out_dir}ld_hapne_pop_trajectory.png",
    script:
        "scripts/fit_hapne_ld.py"


# IBD-based pipeline
rule download_hapibd:
    output:
        "resources/hap-ibd.jar",
    shell:
        """
        curl -L -o {output} https://faculty.washington.edu/browning/hap-ibd.jar
        """


rule download_merge_ibd:
    output:
        "resources/merge-ibd.jar",
    shell:
        """
        curl -L -o {output} https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar
        """


rule run_hapibd:
    input:
        vcf=f"{data_dir}{{name}}.vcf.gz",
        map=f"{data_dir}{{name}}{suffix_map}",
        jar="resources/hap-ibd.jar",
    output:
        "steps/{name}.ibd.gz",
        "steps/{name}.hbd.gz",
    log:
        "logs/hapibd/{name}.log",
    shadow:
        "minimal"
    threads: 2
    conda:
        "conda_environment.yaml"
    params:
        params=config.get("params", ""),
    shell:
        """
        # java -jar hap-ibd.jar gt=$file map=plink.chr$CHR.GRCh38.map  out=IBD/$PREFIX
        java -jar {input.jar} \
            gt={input.vcf} map={input.map} out=steps/{wildcards.name} \
            nthreads={threads} {params} &> /dev/null
        mv steps/{wildcards.name}.log {log}
        """


# According to the HapNe documentation, they recomment to merge ibd and hbd files
# and merge adjacent segments.
rule post_processing_ibd:
    input:
        vcf=f"{data_dir}{{name}}.vcf.gz",
        map=f"{data_dir}{{name}}{suffix_map}",
        ibd="steps/{name}.ibd.gz",
        hbd="steps/{name}.hbd.gz",
        merge_jar="resources/merge-ibd.jar",
    output:
        "steps/{name}.postprocessed.ibd.gz",
    log:
        "logs/postprocessing_ibd/{name}.log",
    shadow:
        "minimal"
    params:
        gap=config.get("gap", 0.6),  # in cM
        discord=config.get("discord", 1),  # at most one discordant homozygote
    conda:
        "conda_environment.yaml"
    shell:
        """
        echo "Processing {wildcards.name}" > {log}
        echo "Processing IBD file" >> {log}
        gunzip -c {input.ibd} | tee >(wc -l > {log}.ibd_lines) | \
            java -jar {input.merge_jar} {input.vcf} {input.map} \
            {params.gap} {params.discord} > steps/{wildcards.name}.postprocessed.ibd 2>> {log}
        echo "Processing HBD file" >> {log}
        gunzip -c {input.hbd} | tee >(wc -l > {log}.hbd_lines) | \
            java -jar {input.merge_jar} {input.vcf} {input.map} \
            {params.gap} {params.discord} >> steps/{wildcards.name}.postprocessed.ibd 2>> {log}
        # Count final output lines
        wc -l steps/{wildcards.name}.postprocessed.ibd > {log}.final_lines
        # Log everything
        echo "IBD lines: $(cat {log}.ibd_lines)" >> {log}
        echo "HBD lines: $(cat {log}.hbd_lines)" >> {log}
        echo "Final lines: $(cat {log}.final_lines)" >> {log}
        # Compression
        gzip steps/{wildcards.name}.postprocessed.ibd
        echo "Done" >> {log}
        # Clean up temp files
        rm {log}.ibd_lines {log}.hbd_lines {log}.final_lines
        """


rule build_histogram:
    input:
        "steps/{name}.postprocessed.ibd.gz",
    output:
        "steps/{name}.ibd.hist",
    log:
        "logs/hist_ibd/{name}.log",
    threads: 1
    params:
        # From the source code: column_cm_length is set to 8 because it's the
        # index of the column that contains the length of the IBD files
        column_cm_length=8,
    conda:
        "conda_environment.yaml"
    shell:
        """
        # Adapted from the source code of HapNe
        gunzip -c {input} | \
            awk -F"\t" '{{l=sprintf("%d", 2*$"{params.column_cm_length}"); c[l]++;}} \
            END {{ for (i=1; i<=40; i++) print i/2/100 "\t" (i+1)/2/100 "\t" 0+c[i]; }}' \
            > {output}
        """


rule fit_hapne_ibd:
    input:
        hists=expand("steps/{name}.ibd.hist", name=names),
        genome_build=config.get("genome_build", ""),
        vcfs=expand(f"{data_dir}{{name}}.vcf.gz", name=names),
    output:
        table=f"{out_dir}ibd_hapne_estimate.csv",
        summary=f"{out_dir}ibd_hapne_summary.txt",
        residuals=f"{out_dir}ibd_hapne_residuals.png",
        popsize=f"{out_dir}ibd_hapne_pop_trajectory.png",
    log:
        "logs/hapne_ibd.log",
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
        nb_bootstraps=config.get("nb_bootstraps"),
        model=config.get("mode"),
    conda:
        "conda_environment.yaml"
    script:
        "scripts/fit_hapne_ibd.py"
