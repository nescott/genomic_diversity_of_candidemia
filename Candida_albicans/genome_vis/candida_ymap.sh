#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:30:00
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --array=2-100
set -ue
set -o pipefail

module load samtools/1.10
module load R/4.1.0

line=${SLURM_ARRAY_TASK_ID}
ref=sc5314  # short ID
snp_bam_file=bam.files # text list of bam files with path
depth_bam_file=gc_corrected.txt # text list of GC-corrected bams with path
fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta

# Get rid of big intermediate files
function finish {
  if [ -e alleles/"${snp_strain}".pileup ]; then
    rm alleles/"${snp_strain}".pileup
  fi
}
trap finish EXIT

# Check for/create output directories
arr=("alleles" "depth" "plots")
for d in "${arr[@]}"; do
    mkdir -p "$d"
done

# Depth and allele counts per bam file from samtools

snp_bam=$(awk -v val="$line" 'NR == val {print $0}' $snp_bam_file)
depth_bam=$(awk -v val="$line" 'NR == val {print $0}' $depth_bam_file)
snp_strain=$(basename "$snp_bam" | cut -d "_" -f 1)
if [ "$snp_strain" == "AMS" ]; then
    snp_strain=$(basename "$snp_bam" | cut -d "_" -f 1,2)
fi
depth_strain=$(basename "$depth_bam" | cut -d "_" -f 1)
if [ "$depth_strain" == "AMS" ]; then
    depth_strain=$(basename "$depth_bam" | cut -d "_" -f 1,2)
fi

# If the snp strain and depth strain don't match, error out
[[ "$depth_strain" != "$snp_strain" ]] && { >&2 echo "Depth and snp IDs don't match"; exit 1; }

# Read depth per position
samtools depth -aa -o depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt "${depth_bam}"

# Add a standard header
sed -i "1s/^/chr\tpos\tdepth\n/" depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt

# Alleles per position (filters on mapping quality of 60, assuming BWA MEM
# alignment), grab the columns of interest
samtools mpileup -q 60 -f "${fasta}" "${snp_bam}" | awk '{print $1, $2, $3, $4, $5}' > alleles/"${snp_strain}".pileup

# Use a borrowed python script to parse mpileup into neat columns of A, T, G, C
python3 /home/selmecki/shared/software/berman_count_snps_v5.py alleles/"${snp_strain}".pileup > alleles/"${snp_strain}"_"${ref}"_putative_SNPs.txt

# Add standard header
sed -i "1s/^/chr\tpos\tref\tA\tT\tG\tC\n/" alleles/"${snp_strain}"_"${ref}"_putative_SNPs.txt

# OPTIONAL - run R script and output plots automatically
# Make sure you've set all the variables in the  R script first!
Rscript --vanilla genome_vis.R depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt alleles/"${snp_strain}"_"${ref}"_putative_SNPs.txt "${depth_strain}"
