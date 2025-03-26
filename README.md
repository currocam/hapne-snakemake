# Hap-Ne snakemake

[HapNe](https://github.com/PalamaraLab/HapNe) is a tool for haplotype-based inference of recent effective population size. To this day, HapNe is only available for human datasets, as it requires the input data to come from a human reference genome.

This repository contains a snakemake workflow that adapts HapNe to work with non-human datasets. Note that I am not the author of HapNe.

Reference: [Haplotype-based inference of recent effective population size in modern and ancient DNA samples](https://doi.org/10.1038/s41467-023-43522-6)

## TL;DR

You can run the two modes of HapNe (LD and IBD) with a toy example by running the following commands:

```bash
snakemake --use-conda -c4 --configfile example/config_ld.yaml
snakemake --use-conda -c4 --configfile example/config_ibd.yaml
```

## Installation

To install the workflow, clone the repository and install the required packages:

```bash
mamba create -n hapne-snakemake -f conda_environment.yaml
conda activate hapne-snakemake
```

Alternatively, you can rely on the `--use-conda` to automatically create the environment if you have snakemake installed.

## HapNe-LD

### Input data

The input data for HapNe-LD is a series of VCF files and genetic maps in SHAPEIT format. Every pair of VCF and genetic map should have the same prefix and correspond to a chromosome (or chromosome arm). An example of the expected input data is in the `example` directory.

### Configuration

The pipeline is configured through a `config.yaml` file. A minimal configuration file is shown below. Additional options can be added such as `apply_filter` or `nb_bootstraps`. Refer to the HapNe documentation for more information.

```yaml
data_dir: "example/"
out_dir: "results/"
map_file_suffix: ".shapeit.map"
method: "ld"
```

The `example`directory should contain the input data. The suffix of the compressed VCF files should be `.vcf.gz`. The suffix of the genetic map files (in shapeit format) can be customized in the `config.yaml` file. Intermediate files will be stored in the `steps` directory. The output of the pipeline is stored in the `results` directory and consists of (1) the estimated _haploid_ effective population size (`*_hapne_estimate.csv`), (2) a summary of the inference (`*_hapne_summary.csv`), (3) a plot of the residuals (`*_hapne_residuals.png`), and (4) a plot of the estimated effective population size (`*_hapne_pop_trajectory.png`).

### Running the pipeline

To run the pipeline, execute the following command:

```bash
snakemake -c8 --configfile config.yaml
```

## HapNe-IBD

### Input data

The input data for HapNe-LD is a series of _phased_ VCF files and genetic maps in _plink_ format. Every pair of VCF and genetic map should have the same prefix and correspond to a chromosome (or chromosome arm).

Additionally, the pipeline requires a TSV file that defines the regions of all analyzed chromosomes (see `example/regions.tsv`). This file should contain the columns 'CHR' (chromosome identifier), 'FROM_BP' (start of the region), 'TO_BP' (end of the region), 'NAME' (name of the region), and 'LENGTH' (length of the region in cM).

### Configuration

The pipeline is configured through a `config.yaml` file. A minimal configuration file is shown below. Additional options can be added such as `apply_filter` or `nb_bootstraps`. Refer to the HapNe documentation for more information.

```yaml
data_dir: "example/"
out_dir: "results/"
map_file_suffix: ".plink.map"
genome_build: "example/regions.tsv"
method: "ibd"
# IBD-detection parameters
hapibd_flags: ""
# Post-processing parameters
gap: 0.6 # Maximum gap between two IBD segments to merge them (in cM)
discord: 1 # at most one discordant homozygote
```

The `example`directory should contain the input data. The suffix of the compressed VCF files should be `.vcf.gz`. The suffix of the genetic map files (in plink format) can be customized in the `config.yaml` file. Intermediate files will be stored in the `steps` directory. The output of the pipeline is stored in the `results` directory and consists of (1) the estimated _haploid_ effective population size (`*_hapne_estimate.csv`), (2) a summary of the inference (`*_hapne_summary.csv`), (3) a plot of the residuals (`*_hapne_residuals.png`), and (4) a plot of the estimated effective population size (`*_hapne_pop_trajectory.png`).

The pipeline will use `HapIBD` to detect IBD segments and process them according to the suggestions of the HapNe authors. The `hapibd_flags` parameter can be used to pass additional flags to `HapIBD`. The `gap` and `discord` parameters are used to merge IBD segments and remove discordant homozygotes, respectively, using the `merge-ibd-segments` tool.

### Running the pipeline

To run the pipeline, execute the following command:

```bash
snakemake -c8 --configfile config.yaml
```

## Analyses of ancient samples

I haven't looked into this, but feel free to make a pull request if you give it a try! Of course, if you have ancient human samples, you should use the original HapNe tool.
