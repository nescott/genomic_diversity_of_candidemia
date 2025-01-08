#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msilarge,msismall
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

export BCFTOOLS_PLUGINS=/home/selmecki/shared/software/software.install/bcftools/1.17/plugins

chr_dir=chr_vcf # location of vcf files to combine
raw_vcf=Calbicans_MEC_strains.vcf.gz  # include .vcf
gene_file=Calbicans_SC5314_A21_sorted_genes.bed.gz
gene_vcf=Calbicans_MEC_info_strains.vcf.gz
regions=SC5314_A21_regions.bed
sorted_vcf=Calbicans_sorted_MEC.vcf.gz # include .vcf
ref=AMS2401
fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/SC5314_A21/C_albicans_SC5314_version_A21-s02-m09-r08_chromosomes.fasta
bcftools_out=Calbicans_MEC_filtered.vcf.gz  # include .vcf
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar
snpeff_config=/home/selmecki/shared/software/snpEff/snpEff.config
snpeff_db=SC5314_s02m09r08  # must be in snpeff.config file and must be name of directory in snpeff data subdir
annotate_vcf=Calbicans_MEC_bwa_filtered_annotated.vcf # include .vcf
genotype_table=Calbicans_MEC_bwa_genotypes.txt  # tab delimited table for use with R scripts (MCA using FactoMineR)

# Load modules
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib/1.9

# Remove intermediate files.
function finish {
  rm samples.txt
  rm chr_files.txt
  rm "${raw_vcf}"*
  rm "${gene_vcf}"*
  rm "${sorted_vcf}"*
  rm "${bcftools_out}"
 }

 trap finish EXIT

# Concatenate regions.
find "$chr_dir" -name "*.vcf" | sort > chr_files.txt
bcftools concat -f chr_files.txt -Oz -o "${raw_vcf}"
bcftools index "${raw_vcf}"

# Add gene name annotation as info column.
bcftools annotate \
    -a "${gene_file}" \
    -c CHROM,FROM,TO,GENE \
    -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
    -o "${gene_vcf}" \
    "${raw_vcf}"

bcftools index "${gene_vcf}"

# Sort samples alphanumerically and remove vars called in SC5314 (AMS2401).
bcftools query -l "${gene_vcf}" | sort > samples.txt
bcftools view -S samples.txt -R "${regions}" -Ou "${gene_vcf}" \
    | bcftools view -e 'GT[0]="alt"' -Ou \
    | bcftools view -s ^"${ref}" -Ou \
    | bcftools view -c 1 -Oz -o "${sorted_vcf}"

bcftools index "${sorted_vcf}"

# Quality and frequency hard filtering.
bcftools view -i \
    "INFO/MQM >39" "${sorted_vcf}" \
    | bcftools view -i "FORMAT/DP[*] > 9" \
    | bcftools view -e 'INFO/TYPE= "complex"' \
    | bcftools norm -f "${fasta}" -m -indels \
    | bcftools +fill-tags - -- -t "FORMAT/VAF" \
    | bcftools view -i "(GT='het' & FORMAT/VAF > 0.14 & FORMAT/VAF < 0.85) | (GT='AA' & FORMAT/VAF > 0.98)" \
    | bcftools view -i "INFO/SAR>=1 & INFO/SAP>0 & INFO/RPL>1 & INFO/RPR>1" \
    | bcftools view -c 1 \
    -Oz -o "${bcftools_out}"

# Annotate using snpeff with manually built database.
java -Xmx9g -jar "${snpeff}" -c "${snpeff_config}" "${snpeff_db}" \
"${bcftools_out}" >  "${annotate_vcf}"

# Subset to biallelic SNPs, output tab-delimited genotype file,
# use for multiple correspondence analysis
module unload bcftools
module load bcftools/1.9

bcftools view -m2 -M2 -v "${annotate_vcf}" \
    | bcftools view -e 'GT="mis"' \
    | bcftools query -H -f '%CHROM\t%POS[\t%GT]\n' >> "${genotype_table}"

sed -i '1s/\[[0-9]\+\]//g' "${genotype_table}"
sed -i '1s/\:GT//g' "${genotype_table}"
sed -i '1s/^\# //' "${genotype_table}"

# Zip and index output.
bgzip "${annotate_vcf}"
tabix "${annotate_vcf}".gz
