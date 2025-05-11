#!/bin/bash
#SBATCH --job-name=UNIREF100
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=cpus
#SBATCH --cpus-per-task=4
##SBATCH --array=1
##SBATCH --ntasks=1
##SBATCH --nodelist=gpu02

wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz