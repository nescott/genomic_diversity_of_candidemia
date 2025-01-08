#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=8gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 45
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-100

set -ue
set -o pipefail

module load jellyfish
module load R/4.1.0

sample_file=../../Calbicans_sequencing_paths.txt
line=${SLURM_ARRAY_TASK_ID}
kmer=21
size=1G
threads=4
read_len=151

# Get vars
sample=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
read1=$(awk -v val="$line" 'NR == val { print $2}' $sample_file)
read2=$(awk -v val="$line" 'NR == val { print $3}' $sample_file)

mkdir -p "${sample}"

jellyfish count -C -m "${kmer}" -s "${size}" -t "${threads}" \
    -o "${sample}"/reads.jf \
    <(zcat "${read1}") <(zcat "${read2}")

jellyfish histo -t "${threads}" "${sample}"/reads.jf > "${sample}"/reads.histo

Rscript genomescope.R "${sample}"/reads.histo  "${kmer}" "${read_len}" "${sample}"
