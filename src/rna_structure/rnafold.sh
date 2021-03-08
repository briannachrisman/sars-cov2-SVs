#!/bin/bash
#SBATCH --job-name=rnafold
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/briannac/logs/rnafold.out
#SBATCH --error=/scratch/groups/dpwall/personal/briannac/logs/rnafold.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15

# sh file is in $MY_SCRATCH/SARS-CoV-2_Microdeletions/src/rna_structure/rnafold.sh #

#VERSION=gisaid_hcov-19_2020_06_03_22

cd $MY_SCRATCH/SARS-CoV-2_Microdeletions
ml python/3.6.1
SVS_TABLE=results/alignment/gisaid_hcov-19_2020_06_03_22.svs_table.tsv

START=1
STOP=29500
FASTA_FOLD=results/rna_structure/NC_045512.2.fasta

PRE_PARAMS="6310 6996 14 RED omark 11074 11541 14 RED omark 20298 22289 14 RED omark 25417 28090 14 RED omark"
RNAplot --pre "${PRE_PARAMS}" --allowFlipping --output-format=ps < ${FASTA_FOLD}
mv NC_045512.2_ss.ps figs/rna_structure/deletion_hotspots.ps


FASTA_FOLD=results/rna_structure/NC_045512.2.fasta
FASTA_OUT=results/rna_structure/${START}_${STOP}.fasta

ANNOTATIONS=results/rna_structure/${START}_${STOP}.annotations.txt 
CONSTRAINTS=results/rna_structure/${START}_${STOP}.constraints.txt 
python3 src/rna_structure/get_fasta_and_annotate.py $START $STOP ${FASTA_OUT} data/reference_genomes/NC_045512.2.fasta ${SVS_TABLE} $ANNOTATIONS $CONSTRAINTS
PRE_PARAMS=$(cat $ANNOTATIONS)
RNAplot --pre "${PRE_PARAMS}" --allowFlipping --output-format=ps < ${FASTA_FOLD}
mv NC_045512.2_ss.ps figs/rna_structure/deletions_all.ps



while read START STOP; do
    echo $START $STOP
    FASTA_OUT=results/rna_structure/${START}_${STOP}.fasta
    FASTA_FOLD=results/rna_structure/NC_045512.2.fasta
    ANNOTATIONS=results/rna_structure/${START}_${STOP}.annotations.txt 
    CONSTRAINTS=results/rna_structure/${START}_${STOP}.constraints.txt 

    python3 src/rna_structure/get_fasta_and_annotate.py $START $STOP ${FASTA_OUT} data/reference_genomes/NC_045512.2.fasta ${SVS_TABLE} $ANNOTATIONS $CONSTRAINTS

    #RNAfold --jobs=15  --commands=$CONSTRAINTS --verbose < ${FASTA_OUT}  > ${FASTA_FOLD}

#/oak/stanford/groups/dpwall/computeEnvironments/ViennaRNA-2.4.14/src/Utils/b2mt.pl results/rna_structure/NC_045512.2.fasta  > results/rna_structure/NC_045512.2_mountain.txt
#awk '{print $1, $NF}' results/rna_structure/NC_045512.2_mountain.txt

#python3 src/rna_structure/get_annotations.py results/alignment/${VERSION}.svs_table.tsv ${ANNOTATIONS}

    PRE_PARAMS=$(cat $ANNOTATIONS)

    RNAplot --pre "${PRE_PARAMS}" --allowFlipping --output-format=ps < ${FASTA_FOLD}

    mv NC_045512.2_ss.ps figs/rna_structure/${START}_${STOP}.ps

done < results/rna_structure/sv_regions.txt
