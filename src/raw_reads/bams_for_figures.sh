#!/bin/bash
#SBATCH --job-name=raw_reads_NCBI
#SBATCH --array=1-13
#SBATCH -p owners
#SBATCH --output=/scratch/users/briannac/logs/raw_reads_NCBI_%a.out
#SBATCH --error=/scratch/users/briannac/logs/raw_reads_NCBI_%a.err
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=50G
######SBATCH --cpus-per-task=1

# File in /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/src/raw_reads/bams_for_figures.sh #
# There are 1179 runs as of right now.

##SLURM_ARRAY_TASK_ID=13
VERSION_NAME=gisaid_hcov-19_seqs


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


for ACC_NUMBER in  SRR11907560 #SRR11780047 SRR11577867 SRR11597222 SRR11621813 SRR11907574 SRR11857975
do

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
    COV2_DEDUP=${INTERMEDIATE_RESULTS}.dedup.bam
    COV2_SORTED=${INTERMEDIATE_RESULTS}.sorted.tmp.bam

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
    samtools index ${COV2_DEDUP}

    ##########################################
    ############### Clean up #################
    ##########################################
    if [ -f "${COV2_DEDUP}" ]; then

        # Remove temporary directories.
        \rm ${INTERMEDIATE_RESULTS}*tmp*

        # Mark pipeline as complete.
        echo "success"
    fi
    echo "done"
done
