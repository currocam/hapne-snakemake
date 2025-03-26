#!/bin/env bash

# Generate datasets
python .test-workflow/simulate_dataset.py .test-workflow/constant_small.yaml .test-workflow/standard_params.yaml example
#python .test-workflow/simulate_dataset.py .test-workflow/constant_big.yaml .test-workflow/big_sample.yaml example
