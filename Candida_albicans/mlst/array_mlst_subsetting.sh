#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=1gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 15
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-100

set -ue
set -o pipefail

module use ~/modulefiles.local

module load bcftools/1.17
module load samtools/1.14

ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
region_file=SC5314_A21_MLST_regions.txt
sample_file=Calbicans_seqs.txt
line=${SLURM_ARRAY_TASK_ID}
vcf=../Calbicans_filter_annotate.vcf.gz

#Read sample file line corresponding to array task ID and get variables
sample=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)

samtools faidx "${ref_fasta}" -r "${region_file}" \
    | bcftools consensus -s "${sample}" -I "${vcf}" > "${sample}"_mlst.fa
