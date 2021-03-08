#!/bin/bash
#SBATCH --job-name=figs
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/briannac/logs/figs.out
#SBATCH --error=/scratch/groups/dpwall/personal/briannac/logs/figs.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15

# sh file is in $MY_SCRATCH/SARS-CoV-2_Microdeletions/src/figs/figs.sh #

ml python/3.6.1
ml viz graphviz

#####################################
### CHANGE DEPENDING ON DATASET #####
VERSION_NAME=gisaid_hcov-19_2020_06_03_22
######################################

PARENT_PATH=/scratch/groups/dpwall/personal/briannac/SARS-CoV-2_Microdeletions

######### Set up file names #########
cd $PARENT_PATH
DELETIONS_TABLE=results/alignment/${VERSION_NAME}.svs_table.tsv
TREE=results/tree/${VERSION_NAME}.tree.nwk
JSON=results/tree/${VERSION_NAME}.nodes.json
MSA_FILEPATH=figs/msa/${VERSION_NAME}
RAW_READS_FILEPATH=figs/raw_reads/${VERSION_NAME}

mkdir ${MSA_FILEPATH}
mkdir ${RAW_READS_FILEPATH}


# Generate deletions for UCSC Genome Browser.
#python3 -u src/figs/generate_bed.py  ${DELETIONS_TABLE} results/UCSC_genome_browser_annotations/${VERSION_NAME}.bed

# Compute histogram.
#python3 -u src/figs/histogram.py results/tree/${VERSION_NAME}.nodes.json results/alignment/${VERSION_NAME}.svs_counter.npy results/alignment/${VERSION_NAME}.alignment_to_ref_pos_dict.npy figs/${VERSION_NAME}_histogram.png

# Compute MSA and raw read figures.
#python3 -u src/figs/msa.py ${DELETIONS_TABLE} results/alignment/${VERSION_NAME}.realigned.fasta NC_045512.2 results/alignment/${VERSION_NAME}.alignment_to_ref_pos_dict.npy ${MSA_FILEPATH}

python3 -u src/figs/raw_reads.py results/alignment/${VERSION_NAME}.raw_reads.tsv data/reference_genomes/NC_045512.2.fasta ${RAW_READS_FILEPATH}


# Draw Tree for deletions.
#mkdir figs/tree/${VERSION_NAME}
#python3 -u src/figs/tree.py ${TREE} ${JSON} results/alignment/${VERSION_NAME}.svs_counter.npy figs/tree/${VERSION_NAME}