# Protein Sequence Search with UniRef100 and HMMER

This repository contains tools for searching similar protein sequences using the [UniRef100](https://www.uniprot.org/help/uniref) database and [jackhmmer](http://hmmer.org/) (part of the HMMER suite). The methodology is based on the [EVcouplings](https://github.com/debbiemarkslab/EVcouplings) framework developed in the Marks Lab for evolutionary sequence analysis and contact prediction.

This implementation specifically applies the methodology described in section 3.1.1 of the supplementary notes of the [EVE (Evolutionary model of Variant Effects) paper](https://www.nature.com/articles/s41586-021-04043-8#Sec10), which details the process of creating high-quality multiple sequence alignments for deep generative modeling of protein sequences.

## HPC Environment

The scripts in this repository, particularly `download-database.sh`, are designed to run in High-Performance Computing (HPC) environments that use the Environment Modules system. The script includes SLURM job scheduler directives for resource allocation and job management.

```bash
#!/bin/bash
#SBATCH --job-name=UNIREF100
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=cpus
#SBATCH --cpus-per-task=4
```

These directives configure the script to run as a SLURM job with appropriate resources for downloading the large UniRef100 database.

## Setup

### Create the Conda Environment

Create a conda environment using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
```

This will create a new environment named `msa` with all necessary dependencies.

### Activate the Environment

```bash
conda activate msa
```

## Database Preparation

### Download UniRef100 Database

Use the included script to download the UniRef100 database:

```bash
bash download-database.sh
```

The script downloads the UniRef100 database (approximately 114GB) from the UniProt FTP server.

### Extract the Database

After downloading, extract the compressed database file:

```bash
gunzip uniref100.fasta.gz
```

## Creating Multiple Sequence Alignments

This repository provides a SLURM job script `create_msa.sh` that handles the creation of Multiple Sequence Alignments (MSAs) in an HPC environment that uses environment modules. The script automatically loads the necessary modules and runs the protein sequence search against the UniRef100 database.

### Using create_msa.sh

The script is configured to run in an HPC environment with SLURM job scheduler and environment modules:

```bash
#!/bin/bash
#SBATCH --job-name=EVcouplings
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=cpus
#SBATCH --cpus-per-task=32 // using 32 cores


ml load anaconda3/2024.10

conda activate evecp

python main.py --sequence_id Q1KHJ8_9INFA --outdir data/MSA/tQ1KHJ8_9INFA
```

To run the script with your protein of interest:

1. Edit the script to change the `--sequence_id` parameter to your protein's UniProt ID
2. Modify the output directory as needed
3. Submit the job to the SLURM scheduler:

```bash
sbatch create_msa.sh
```

The script will:
1. Load the anaconda3 module in your HPC environment
2. Activate the evecp conda environment
3. Run the MSA creation pipeline for your specified protein

## MSA Creation for EVE Model

This repository provides scripts necessary for creating Multiple Sequence Alignments (MSAs) required by the [EVE (Evolutionary model of Variant Effects) model](https://github.com/OATML-Markslab/EVE). The MSAs generated with these tools serve as input data for EVE, which predicts the pathogenicity of protein variants using deep generative models of evolutionary data.

### EVE MSA Requirements

The EVE model requires MSAs with:
- Sequences that align to at least 50% of the target protein sequence
- Columns with at least 70% residue occupancy
- Sufficient depth of sequences (ideally between 10L and 100,000 sequences, where L is the protein length)

### Creating MSAs for EVE

The `create_msa.sh` script handles the creation of MSAs suitable for EVE analysis through EVcouplings. The script automates the process of:

1. Preparing the query protein sequence
2. Running jackhmmer against UniRef100 with appropriate parameters
3. Processing and filtering the resulting alignments according to EVE requirements

For more details about EVE's approach to MSA creation, see the [EVE GitHub repository](https://github.com/OATML-Markslab/EVE).


## Additional Resources

- [HMMER Documentation](http://hmmer.org/documentation.html)
- [UniProt/UniRef100](https://www.uniprot.org/help/uniref)
- [EVcouplings Documentation](https://evcouplings.org)
- [EVE Model](https://evemodel.org/)

## Note

The UniRef100 database is very large (approximately 114GB compressed). Ensure you have sufficient disk space before downloading. 