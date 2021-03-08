#!/bin/bash
#SBATCH --job-name=align
#SBATCH -p dpwall
#SBATCH --output=/scratch/users/briannac/logs/align.out
#SBATCH --error=/scratch/users/briannac/logs/align.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=120G
#SBATCH --cpus-per-task=15

# sh file is in $MY_HOME/SARS-CoV-2_structural_variations/src/alignment/align.sh #


ml python/3.6.1


PARENT_PATH=$MY_HOME/SARS-CoV-2_structural_variations

######### Set up file names #########
cd $PARENT_PATH

#####################################
### CHANGE DEPENDING ON DATASET #####
VERSION_NAME=gisaid_seqs
NCBI_READS_TABLE=data/raw_reads/ncbi_sra_metadata.txt
######################################



RESULTS=results/alignment
INTERMEDIATE_FILES=intermediate_files/alignment
REF=data/reference_genome/NC_045512.2.fasta
SEQS_RAW=$INTERMEDIATE_FILES/${VERSION_NAME}.fasta
SEQS_FILT=$INTERMEDIATE_FILES/${VERSION_NAME}.filt.fasta
SEQS_W_REF=$INTERMEDIATE_FILES/${VERSION_NAME}.w_ref.fasta
SEQS_NUM=$INTERMEDIATE_FILES/${VERSION_NAME}.num.fasta
ALIGNMENT=$INTERMEDIATE_FILES/${VERSION_NAME}.fasta

REALIGNMENT=$RESULTS/${VERSION_NAME}.realigned.fasta
SV_COUNTER=$RESULTS/${VERSION_NAME}.svs_counter.npy
SV_DICT=$RESULTS/${VERSION_NAME}.svs_dict.npy
SV_TABLE=$RESULTS/${VERSION_NAME}.svs_table.tsv
FILT_READS_TABLE=$RESULTS/${VERSION_NAME}.raw_reads.tsv
REF_DICT=$RESULTS/${VERSION_NAME}.alignment_to_ref_pos_dict.npy
#SV_FINAL_COUNTER=$RESULTS/${VERSION_NAME}.svs_counter_final.npy
#SV_FINAL_DICT=$RESULTS/${VERSION_NAME}.svs_dict_final.npy
#SV_FINAL_TABLE=$RESULTS/${VERSION_NAME}.svs_table_final.tsv

# Filter sequences.
echo "Filtering...."

cat intermediate_files/alignment/gisaid_hcov* > ${SEQS_RAW}

python3 src/alignment/filter_fasta.py ${SEQS_RAW} ${SEQS_FILT}

cat ${REF} ${SEQS_FILT} > ${SEQS_W_REF}

awk '{if (/^>/) print ">seq_"(++i)"|" substr($0,2); else print $0;}' $SEQS_W_REF > $SEQS_NUM

echo "Aligning...."
augur align \
  --sequences $SEQS_NUM \
  --output $ALIGNMENT \
  --nthreads=auto
  
echo "Finished Aligning"

# Counts & realigns the SVs.
python3 src/alignment/count_svs_and_realign.py $ALIGNMENT ${REALIGNMENT} ${SV_COUNTER} ${SV_DICT} ${REF_DICT} 'NC_045512.2'

# Computes the IDs of the SRI experiments with indels, creates a table.
python3 src/alignment/raw_reads_with_svs.py ${SV_DICT} ${SV_TABLE} ${NCBI_READS_TABLE} ${FILT_READS_TABLE}
