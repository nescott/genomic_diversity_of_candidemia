#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 50
#SBATCH -p msilarge,msismall
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

out_file=Calbicans_MEC_bwa_filtered_annotated

# load modules, use local for newer bcftools
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib/1.9
module load gatk/4.1.2

export BCFTOOLS_PLUGINS=/home/selmecki/shared/software/software.install/bcftools/1.17/plugins

# remove int files
function finish {
  rm -f "${out_file}".table
  rm -f "${out_file}".ann.table
  rm -f "${out_file}".first_ann.table
  rm -f "${out_file}".tsv
}

trap finish EXIT

# get info in tabular form
gatk VariantsToTable \
    -V "${out_file}".vcf.gz \
    -F CHROM -F POS -F GENE -F REF -F ALT -GF GT -GF VAF \
    -O "${out_file}".table

# grab annotation
bcftools query -f '%INFO/ANN\n' \
    "${out_file}".vcf.gz > "${out_file}".ann.table

# keep only the first
awk 'BEGIN{FS="|"; OFS="\t"} {print $2, $3, $4, $10, $11}' \
    "${out_file}".ann.table > "${out_file}".first_ann.table

# fix the header
sed -i '1i Annotation\tImpact\tGene_Name\tCoding_Change\tAA_change' \
    "${out_file}".first_ann.table

# do row numbers still match?
[[ $(wc -l < $out_file.table) -ne $(wc -l < $out_file.first_ann.table) ]] && { >&2 echo "table lengths don't match"; exit 1; }

# final tsv output
paste "${out_file}".table "${out_file}".first_ann.table > "${out_file}".tsv
tr "\t" "," < "${out_file}".tsv > "${out_file}".csv
