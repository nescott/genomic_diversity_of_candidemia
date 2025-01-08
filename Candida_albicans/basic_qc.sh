#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=4:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail
unset DISPLAY

fastq_dir=/scratch.global/scot0854/calbicans/trimmed_fastq/
bam_dir=bam
temp_dir=/scratch.global/scot0854/calbicans  #scratch directory
file_name="$(date +%F)_Calbicans_qc"

# Use local modules
module use /home/selmecki/shared/software/modulefiles.local

# Load modules
module load fastqc
module load qualimap
module load multiqc

# fastqc on bbduk-trimmed files
find "$fastq_dir" -name "*trimmed_*P.fq" -exec fastqc  -t 8 -o "$temp_dir" {} \;

# bam qc
unset DISPLAY # qualimap won't work on cluster without this
find "$bam_dir" -name "*.bam" \
  -exec qualimap bamqc -bam {} -outdir $temp_dir/{} --java-mem-size=4G  \;

# multiqc
multiqc /scratch.global/scot0854/calbicans logs -n "$file_name"
