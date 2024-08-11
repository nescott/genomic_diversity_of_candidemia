#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=48:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1,3-9

# call variants for all samples in a population using freebayes, subsetting by region
set -ue
set -o pipefail

species=Calbicans
ref=sc5314_a21
ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
line=${SLURM_ARRAY_TASK_ID}
bam_list=bam.files
region_list=Calbicans_chroms.txt

chr=$(awk -v val="$line" 'NR == val { print $0}' $region_list)

#Load modules
module load samtools/1.10
module load freebayes/20180409

freebayes -f "${ref_fasta}" \
  -C 5 \
  -p 2 \
  -r "${chr}" \
  -L "${bam_list}" \
  -v "chr_vcf/${species}_${ref}_${chr}.vcf"
