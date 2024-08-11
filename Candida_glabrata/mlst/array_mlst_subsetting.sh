#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 30
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-98

set -ue
set -o pipefail

module use ~/modulefiles.local

module load bcftools/1.17
module load samtools/1.14

ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/Cglabrata/CBS138/GCF_000002545.3_ASM254v2_genomic.fna.gz
region_file=GCF_000002545.3_samtools_regions.txt
sample_file=../../cglab_todo.txt
line=${SLURM_ARRAY_TASK_ID}
vcf=../Cglabrata_filter_annotate.vcf.gz

#Read sample file line corresponding to array task ID and get variables
sample=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)

samtools faidx "${ref_fasta}" -r "${region_file}" \
    | bcftools consensus -s "${sample}" -I "${vcf}" > "${sample}"_mlst.fa
