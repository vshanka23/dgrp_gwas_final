#!/bin/bash
#
#SBATCH --job-name=mafft_mito
#SBATCH --cpus-per-task=40
#SBATCH --partition=bigmem
#SBATCH --time=30-00:00:00
#SBATCH --mem=700g
#SBATCH --output=/scratch2/vshanka/Projects/vshanka_dgrp3_dna/8_mt_hap/mt_pipeline_results/logs/mafft_output.%j.out
#SBATCH --error=/scratch2/vshanka/Projects/vshanka_dgrp3_dna/8_mt_hap/mt_pipeline_results/logs/mafft_error.%j.err
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

ml mafft

cd /scratch2/vshanka/Projects/vshanka_dgrp3_dna/8_mt_hap/mt_pipeline_results/

mafft --auto --thread 40 mt_consensus_all.fasta > mt_consensus_all_aligned.fasta