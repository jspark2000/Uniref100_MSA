#!/bin/bash
#SBATCH --job-name=EVcouplings
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=cpus
#SBATCH --cpus-per-task=32

ml load anaconda3/2024.10

conda activate evecp

python main.py --sequence_id Q1KHJ8_9INFA --outdir data/MSA/tQ1KHJ8_9INFA