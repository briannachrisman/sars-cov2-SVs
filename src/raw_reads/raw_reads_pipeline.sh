#!/bin/bash
#SBATCH --job-name=raw_reads_NCBI
#SBATCH --array=1
#SBATCH -p owners
#SBATCH --output=/scratch/users/briannac/logs/raw_reads_NCBI_%a.out
#SBATCH --error=/scratch/users/briannac/logs/raw_reads_NCBI_%a.err
#SBATCH --time=00:40:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=100G

# File in /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/src/raw_reads/raw_reads_pipeline.sh #
# There are 1179 runs as of right now.

##SLURM_ARRAY_TASK_ID=13
VERSION_NAME=gisaid_hcov-19_seqs
RAW_READ_METADATA=/home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/intermediate_files/raw_reads/ncbi_metadata_not_finished.tsv
#results/alignment/${VERSION_NAME}.raw_reads.tsv

########## COMMENT OUT AFTER TESTING!!!!! ###############
N=$((SLURM_ARRAY_TASK_ID + 0)) 

# Module load.
ml biology
ml samtools
module load gatk/4.1.4.1

# Accesion ID list downloaded by from NCBI SRA reader.

##########################################
############ Set up read data. ###########
##########################################

# Set variable names.
PARENT_DIR=$MY_HOME/SARS-CoV-2_structural_variations
cd $PARENT_DIR

#data/reference_genomes/NC_045512.2
COV2_REF=data/reference_genome/NC_045512.2.fasta

head -${N} ${RAW_READ_METADATA} | tail -n 1 > $MY_SCRATCH/tmp/test_${N}.tmp

while read -r ACC_NUMBER_ ASSAY_ PLATFORM_ PAIRED_; do
    ACC_NUMBER=$ACC_NUMBER_
    ASSAY=$ASSAY_
    PLATFORM=$PLATFORM_
done < $MY_SCRATCH/tmp/test_${N}.tmp
\rm $MY_SCRATCH/tmp/test_${N}.tmp

echo "###############################"
echo "ACC_NUMBER: " $ACC_NUMBER
echo "ASSAY: " $ASSAY
echo "PLATFORM: " $PLATFORM
echo "###############################"


##########################################
######## SETUP FILE DIRECTORIES ##########
##########################################
INTERMEDIATE_RESULTS=/home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/intermediate_files/raw_reads/${ACC_NUMBER}
FINAL_RESULTS=results/raw_reads/${ACC_NUMBER}
R1=${INTERMEDIATE_RESULTS}_1.tmp.fastq
R=${INTERMEDIATE_RESULTS}.tmp.fastq
R2=${INTERMEDIATE_RESULTS}_2.tmp.fastq
R1_clean=${INTERMEDIATE_RESULTS}_1_clean.tmp.fastq
R2_clean=${INTERMEDIATE_RESULTS}_2_clean.tmp.fastq
R_clean=${INTERMEDIATE_RESULTS}_clean.tmp.fastq

COV2_ALL=${INTERMEDIATE_RESULTS}.cov2.tmp.bam
COV2_UNMAPPED=${INTERMEDIATE_RESULTS}.unmapped.tmp.bam
COV2_MAPPED=${INTERMEDIATE_RESULTS}.mapped.tmp.bam
COV2_GOOD_MAP=${INTERMEDIATE_RESULTS}.good.tmp.bam
COV2_DEDUP=${INTERMEDIATE_RESULTS}.markdup.tmp.bam
COV2_DEDUP_METRICS=${INTERMEDIATE_RESULTS}.tmp.dupmetrics.tmp.txt
COV2_INDEL=${INTERMEDIATE_RESULTS}.indel.tmp.bam
COV2_SORTED=${INTERMEDIATE_RESULTS}.sorted.tmp.bam
COV2_VITERBI=${INTERMEDIATE_RESULTS}.viterbi.tmp.bam
COV2_VITERBI_SORTED=${INTERMEDIATE_RESULTS}.tmp.bam
COV2_VCF=${INTERMEDIATE_RESULTS}.tmp.cov.vcf
VCF_FINAL=${INTERMEDIATE_RESULTS}.cov.vcf

##########################################
########### Download read data. ##########
##########################################
echo "Downloading data...."

############ SRA READS ##################
fastq-dump --split-3  --skip-technical ${ACC_NUMBER}
mv ${ACC_NUMBER}.fastq ${R}
PAIRED=SINGLE
if [ -f ${ACC_NUMBER}_2.fastq ]; then
    PAIRED=PAIRED
    mv ${ACC_NUMBER}_1.fastq ${R1}
    mv ${ACC_NUMBER}_2.fastq ${R2}
fi
echo "downloaded data"
echo ${PAIRED}

##########################################
############# QUALITY CONTROL ############
##########################################
echo "quality control...."
if [ ${PAIRED} = "SINGLE" ]; then
    fastp -i ${R} \
    -o ${R_clean} \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 20 \
    --length_required 50 
fi

if [ ${PAIRED} = "PAIRED" ]; then
    fastp -i ${R1} -I ${R2} \
    -o ${R1_clean} -O ${R2_clean} \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 20 \
    --length_required 50 
fi

echo "done with quality control"


##########################################
############ Align to SARS-CoV-2. ########
##########################################
# Index fastq file (this takes a while.)
# Align to SARS-CoV-2 reference.
# Make files of mapped and unmapped reads.
# Sort mapped reads for pileup later on.
echo "aligning to sars-cov-2..."
if [ ${PAIRED} = "SINGLE" ]; then
    bwa mem ${COV2_REF} ${R_clean} -o ${COV2_ALL} -t 5
fi
if [ ${PAIRED} = "PAIRED" ]; then
    bwa mem ${COV2_REF} ${R1_clean} ${R2_clean} -o ${COV2_ALL} -t 5
fi
echo "aligned to sars-cov-2."

samtools view  -f 0x04 -b -o ${COV2_UNMAPPED} ${COV2_ALL}

echo "# unaligned to cov2: " $(samtools view -c ${COV2_UNMAPPED})

if [ ${PAIRED} = "SINGLE" ]; then
    samtools view  -F 0x04 -b -o ${COV2_MAPPED} ${COV2_ALL}
fi

if [ ${PAIRED} = "PAIRED" ]; then
    samtools view  -f 0x02 -b -o ${COV2_MAPPED} ${COV2_ALL}
fi

echo "# aligned to cov2: " $(samtools view -c ${COV2_MAPPED})

samtools sort -o ${COV2_SORTED} ${COV2_MAPPED} 
samtools index ${COV2_SORTED}

##########################################
#### Mark and delete PCR repilcates ######
##########################################
echo "calling variants..."
gatk MarkDuplicates -I ${COV2_SORTED} -O ${COV2_DEDUP} -M ${COV2_DEDUP_METRICS} --REMOVE_DUPLICATES

lofreq indelqual --ref ${COV2_REF} -o ${COV2_INDEL} --dindel --verbose  ${COV2_DEDUP} 
lofreq viterbi --ref ${COV2_REF} -o ${COV2_VITERBI} --verbose ${COV2_INDEL} 

samtools sort -o ${COV2_VITERBI_SORTED} ${COV2_VITERBI} 
samtools index ${COV2_VITERBI_SORTED}

lofreq call --ref ${COV2_REF} -o ${COV2_VCF} --call-indels --verbose ${COV2_VITERBI_SORTED} 
sed -i "s/\t\.\t/\t${ACC_NUMBER}\t/g" ${COV2_VCF} # Replace '.' id with accession number.

LINES=$(wc -l COV2_VCF | awk '{print $1 }')
if [ $LINES -eq 0 ]
then 
    echo "FAIL"
else
    grep 'NC' ${COV2_VCF}  | tail -n +3 > ${VCF_FINAL}
fi

##########################################
############### Clean up #################
##########################################
if [ -f "${VCF_FINAL}" ]; then

    # Remove temporary directories.
    \rm ${INTERMEDIATE_RESULTS}*tmp*
    
    # Mark pipeline as complete.
    echo "success"
fi
echo "done"



## After done, run:
## cat /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/intermediate_files/raw_reads/*.cov.vcf > /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/results/raw_reads/all.cov.vcf