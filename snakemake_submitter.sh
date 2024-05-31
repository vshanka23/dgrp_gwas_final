#!/bin/bash
#
#SBATCH --job-name=<job_name>
#SBATCH --ntasks=1   
#SBATCH --partition=<controller partition name>
#SBATCH --time=30-00:00:00
#SBATCH --mem=2gb
#SBATCH --output=<path to working directory>/log/test_output_%j.txt
#SBATCH --error=<path to working directory>/log/test_error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=<uid@domain.edu>

cd <path to working directory>/
#mkdir -p ./{log,logs_slurm}

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh
conda activate snakemake

#--dag | display | dot
#-p -n \

#GEMMA with GEMMA GRM
snakemake \
-s Snakefile \
--profile slurm \
--configfile config.yaml \
--latency-wait 120 \
-p

#GEMMA with GCTA GRM (comment out the previous snakemake and uncomment next line if you want to run with GCTA GRM and --make-grm-inbred flag)
#snakemake -s Snakefile_GCTA --profile slurm --configfile config.yaml --latency-wait 120