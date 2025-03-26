#!/usr/bin/env python3
import msprime, demes
import pandas as pd
import os, subprocess, yaml, sys


def main(demes_file: str, config: str, outdir: str):
    # Load the configuration file
    with open(config, "r") as ymlfile:
        cfg = yaml.safe_load(ymlfile)
    num_individuals = cfg["num_individuals"]
    sequence_length = float(cfg["sequence_length"])
    mutation_rate = float(cfg["mutation_rate"])
    recombination_rate = float(cfg["recombination_rate"])
    num_chromosomes = cfg["num_chromosomes"]
    seed = cfg["seed"]
    chrom_names = [str(i) for i in range(num_chromosomes)]
    graph = demes.load(demes_file)
    demography = msprime.Demography.from_demes(graph)
    # Define out files
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    seeds = [seed + i for i in range(num_chromosomes)]
    # Grow region dataframe row by row
    # CHR	FROM_BP	TO_BP	NAME	LENGTH
    regions = []
    for index, seed in enumerate(seeds):
        # Add some noise to the chromosome length
        length = sequence_length + seed * 1e3
        ts = msprime.sim_ancestry(
            samples=num_individuals,
            demography=demography,
            random_seed=seed,
            recombination_rate=recombination_rate,
            sequence_length=length,
        )
        length_in_cm = length * recombination_rate * 100
        regions.append(
            {
                "CHR": chrom_names[index],
                "FROM_BP": 0,
                "TO_BP": length,
                "NAME": chrom_names[index],
                "LENGTH": length_in_cm,
            }
        )
        # Write genetic map file in SHAPEIT format
        mapfile = os.path.join(outdir, f"{chrom_names[index]}.shapeit.map")
        with open(mapfile, "w") as gm:
            gm.write("position\tCOMBINED_rate(cM/Mb)\tGenetic_Map(cM)\n")
            gm.write("0\t1.0\t0\n")
            gm.write(f"{int(length)}\t1.0\t{length_in_cm}\n")
        # Write genetic map file in PLINK format
        mapfile = os.path.join(outdir, f"{chrom_names[index]}.plink.map")
        with open(mapfile, "w") as gm:
            # 1  rs100  0  1\n1  rs101  1  1000000
            gm.write(f"{chrom_names[index]} rs100 0 1\n")
            gm.write(f"{chrom_names[index]} rs101 {length_in_cm} {int(length)}\n")

        mts = msprime.sim_mutations(ts, rate=mutation_rate)
        n_dip_indv = int(mts.num_samples / 2)
        indv_names = [f"tsk_{i}indv" for i in range(n_dip_indv)]
        outfile = os.path.join(outdir, f"{chrom_names[index]}.vcf")
        with open(outfile, "w") as vcf_file:
            mts.write_vcf(vcf_file, individual_names=indv_names)
        with open("chr_name_conv.txt", "a") as chr_conv:
            chr_conv.write(f"1 {chrom_names[index]}\n")
        result = subprocess.run(
            f"bcftools annotate --rename-chrs chr_name_conv.txt {outfile} | bgzip > {outfile}.gz",
            shell=True,
        )
        if result.returncode != 0:
            print(f"Error running bcftools for {chrom_names[index]}")
        os.remove(outfile)
        os.remove("chr_name_conv.txt")
    # Write region file
    region_df = pd.DataFrame.from_dict(regions)
    region_df.to_csv(os.path.join(outdir, "regions.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) == 4:
        demes_file = sys.argv[1]
        config = sys.argv[2]
        outdir = sys.argv[3]
        main(demes_file, config, outdir)
