#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=3:00:00
#SBATCH -p amdsmall,amdlarge,amd512
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-104

set -ue
set -o pipefail

species=Cglabrata
ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/Cglabrata/CBS138_s05-m03-r02/C_glabrata_CBS138_version_s05-m03-r02_chromosomes.fasta.gz
sample_file=Cglabrata_seq_paths.txt
line=${SLURM_ARRAY_TASK_ID}

#Read sample file line corresponding to array task ID and get variables
read1=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
read2=$(awk -v val="$line" 'NR == val { print $2}' $sample_file)
strain=$(basename $read1 | cut -d "_" -f 1)

# Load modules for trimming and aligning
module use /home/selmecki/shared/software/modulefiles.local

module load bwa/0.7.17
module load samtools/1.10
module load bbmap

mkdir -p trimmed_fastq logs bam

# Adapter and quality trimming using JGI BBTools data preprocessing guidelines
## trim adapters
bbduk.sh in1="${read1}" in2="${read2}" \
out1=trimmed_fastq/"${strain}"_trim_adapt1.fq out2=trimmed_fastq/"${strain}"_trim_adapt2.fq \
ref=adapters ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo

## contaminant (phix) filtering
bbduk.sh in1=trimmed_fastq/"${strain}"_trim_adapt1.fq in2=trimmed_fastq/"${strain}"_trim_adapt2.fq \
out1=trimmed_fastq/"${strain}"_unmatched1.fq out2=trimmed_fastq/"${strain}"_unmatched2.fq \
outm1=trimmed_fastq/"${strain}"_matched1.fq outm2=trimmed_fastq/"${strain}"_matched2.fq \
ref=phix,artifacts k=31 hdist=1 stats=logs/"${strain}"_phistats.txt

## quality trimming (bbduk user guide recommends this as separate step from adapter trimming)
bbduk.sh in1=trimmed_fastq/"${strain}"_unmatched1.fq in2=trimmed_fastq/"${strain}"_unmatched2.fq \
out1=trimmed_fastq/"${strain}"_trimmed_1P.fq out2=trimmed_fastq/"${strain}"_trimmed_2P.fq qtrim=rl trimq=10

# Reference alignment, fix mate-pair errors from alignment, sort, mark duplicates
# Including sample information to ensure unique read groups if freebayes is used
bwa mem -t 8 -R "@RG\tID:${species}_${strain}\tPL:ILLUMINA\tSM:${strain}" \
"${ref_fasta}" trimmed_fastq/"${strain}"_trim_sort_1P.fq \
trimmed_fastq/"${strain}"_trim_sort_2P.fq \
| samtools fixmate -m - - \
| samtools sort -l 0 -T "${strain}" -@8 - \
| samtools markdup -@8 - bam/"${strain}"_trimmed_bwa_sorted_markdup.bam

#reindex
samtools index bam/"${strain}"_trimmed_bwa_sorted_markdup.bam

#basic stats
samtools flagstat bam/"${strain}"_trimmed_bwa_sorted_markdup.bam \
> logs/"${strain}"_trimmed_bwa_sorted_markdup.stdout
