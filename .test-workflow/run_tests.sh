#!/bin/env bash
# Fail on error
set -e
echo "Running LD tests"
snakemake -c4 -F --configfile example/config_ld.yaml
snakemake -c4 -F --configfile .test-workflow/config_ld2.yaml
echo "Running IBD tests"
snakemake -c4 -F --configfile example/config_ibd.yaml
snakemake -c4 -F --configfile .test-workflow/config_ibd2.yaml
