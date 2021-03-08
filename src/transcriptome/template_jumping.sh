#!/bin/bash
#SBATCH --job-name=template_jumping
#SBATCH -p dpwall
#SBATCH --output=/scratch/users/briannac/logs/template_jumping.out
#SBATCH --error=/scratch/users/briannac/logs/template_jumping.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15

# sh file is in $MY_HOME/SARS-CoV-2_structural_variations/src/transcriptome/template_jumping.sh #

ml python/3.6.1
ml py-pysam/0.15.3_py36

PARENT_PATH=/home/groups/dpwall/briannac/SARS-CoV-2_structural_variations
cd $PARENT_PATH

# Download transcriptome dataset.
#wget https://osf.io/2wktx/download 
#mv download intermediate_files/transcriptome/VeroInf24h.viral_genome.bam

#samtools sort -o intermediate_files/transcriptome/VeroInf24h.viral_genome_sorted.bam #intermediate_files/transcriptome/VeroInf24h.viral_genome.bam 
#samtools index intermediate_files/transcriptome/VeroInf24h.viral_genome_sorted.bam
######### Set up file names #########

python3 -u src/transcriptome/template_jumping.py intermediate_files/transcriptome/VeroInf24h.viral_genome_sorted.bam results/transcriptome/read_ends.csv