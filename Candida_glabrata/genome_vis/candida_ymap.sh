#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:30:00
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=38
set -ue
set -o pipefail

module load samtools/1.10
module load R/4.1.0

line=${SLURM_ARRAY_TASK_ID}
ref=cgd_s05m03r02
depth_bam_file=gc_corrected.txt # text list of GC-corrected bams with path
feature_file=Cglabrata_cgd_s05m03r02_features.txt
chr_labels=Cglabrata_chr_labels.txt

mkdir -p depth

depth_bam=$(awk -v val="$line" 'NR == val {print $0}' $depth_bam_file)
depth_strain=$(basename "$depth_bam" | cut -d "_" -f 1)
if [ "$depth_strain" == "AMS" ]; then
    depth_strain=$(basename "$depth_bam" | cut -d "_" -f 1,2)
fi

# Read depth per position
samtools depth -aa -o depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt "${depth_bam}"

# Add a standard header
sed -i "1s/^/chr\tpos\tdepth\n/" depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt

# Run R script and output plots automatically
Rscript --vanilla genome_vis.R \
    depth/"${depth_strain}"_"${ref}"_gc_corrected_depth.txt \
    "${depth_strain}" \
    "${feature_file}" \
    "${chr_labels}"
