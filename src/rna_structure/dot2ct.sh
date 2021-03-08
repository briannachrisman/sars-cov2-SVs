#!/bin/bash
#SBATCH --job-name=dot2ct
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/briannac/logs/dot2ct.out
#SBATCH --error=/scratch/groups/dpwall/personal/briannac/logs/dot2ct.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1

# sh file is in $MY_SCRATCH/SARS-CoV-2_Microdeletions/src/rna_structure/dot2ct.sh #

dot2ct $MY_SCRATCH/SARS-CoV-2_Microdeletions/results/rna_structure/NC_045512.2.fasta $MY_SCRATCH/SARS-CoV-2_Microdeletions/results/rna_structure/NC_045512.2.ct
