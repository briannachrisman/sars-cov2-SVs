#!/bin/bash
#SBATCH --job-name=spot_rna
#SBATCH -p bigmem
#SBATCH --output=/scratch/users/briannac/logs/spot_rna.out
#SBATCH --error=/scratch/users/briannac/logs/spot_rna.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=2T

# sh file is in $MY_HOME/SARS-CoV-2_structural_variations/src/rna_structure/spot_rna.sh #

python3.6 /oak/stanford/groups/dpwall/computeEnvironments/SPOT-RNA/SPOT-RNA.py  \
--inputs /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/data/reference_genome/NC_045512.2.fasta \
--outputs /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/results/rna_structure/spot_rna