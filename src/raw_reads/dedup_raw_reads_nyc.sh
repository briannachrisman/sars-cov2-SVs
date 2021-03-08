#!/bin/bash
## Removing duplicates from NYC data & copying to interesting directory.

IDX=V20_0333
samtools view -H /scratch/groups/dpwall/personal/briannac/covid19_evolution/data/raw_fastqs/raw_reads_private/SARS_CoV2_RawReads_MauranoMatt/${IDX}*.bam | sed -e 's/NC_045512v2/NC_045512.2/'| samtools reheader - /scratch/groups/dpwall/personal/briannac/covid19_evolution/data/raw_fastqs/raw_reads_private/SARS_CoV2_RawReads_MauranoMatt/${IDX}*.bam \
> /scratch/groups/dpwall/personal/briannac/covid19_evolution/results/raw_reads_with_deletion/${IDX}.tmp.bam
samtools index /scratch/groups/dpwall/personal/briannac/covid19_evolution/results/raw_reads_with_deletion/${IDX}.tmp.bam
cd /scratch/groups/dpwall/personal/briannac/covid19_evolution/results/raw_reads_with_deletion

module load gatk
#samtools sort -o ${IDX}.sorted.tmp.bam ${IDX}.tmp.bam 
gatk MarkDuplicates -I ${IDX}.tmp.bam -O ${IDX}.dedup.bam -M ${IDX}.tmp.dedupmetrics --REMOVE_DUPLICATES

samtools index ${IDX}.dedup.bam
\rm *.tmp*


