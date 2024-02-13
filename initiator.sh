#!/bin/bash
mkdir -p ./{log,logs_slurm} | sbatch ./snakemake_submitter.sh
