#!/bin/bash
#SBATCH --job-name=rna_mxfold
#SBATCH -p dpwall
#SBATCH --output=/scratch/users/briannac/logs/rna_mxfold.out
#SBATCH --error=/scratch/users/briannac/logs/rna_mxfold.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=50G

# sh file is in $MY_HOME/SARS-CoV-2_structural_variations/src/rna_structure/mxfold.sh #

/oak/stanford/groups/dpwall/computeEnvironments/mxfold/mxfold/build/mxfold /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/data/reference_genome/NC_045512.2.fasta  \
>  /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/results/rna_structure/mxfold.ct 