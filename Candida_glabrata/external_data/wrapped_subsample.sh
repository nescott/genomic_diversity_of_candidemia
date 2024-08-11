#!/bin/bash
set -ue
set -o pipefail

line=$1
sample_file=bams_to_subsample.txt

# Read sample file line corresponding to array task ID and get variables
read1=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
strain=$(basename "$read1" | cut -d "_" -f 1)

# Load modules for trimming and aligning
module use /home/selmecki/shared/software/modulefiles.local

module load bbmap
module load samtools/1.10

reformat.sh in="${read1}" \
out=combined_data/subsample_bam/"$strain"_subset.bam \
sampleseed=2 samplereadstarget=5000000

# Reindex
samtools index combined_data/subsample_bam/"${strain}"_subset.bam

# Basic qc
samtools flagstat combined_data/subsample_bam/"${strain}"_subset.bam \
   >  combined_data/logs/"${strain}".log

module load qualimap
unset DISPLAY # Qualimap won't work on cluster without this
qualimap bamqc -bam combined_data/subsample_bam/"${strain}"_subset.bam \
    -outdir combined_data/subsample_bam/"${strain}" --java-mem-size=4G
