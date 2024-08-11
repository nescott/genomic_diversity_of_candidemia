#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 20
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

module load bedtools/2.29.2

breakpoint_file=cnv_breakpoints.bed
repeat_file=selmecki_cgd_unique_repeats.gff3
out_file=Calb_MEC_cnv_closest.txt

bedtools closest -D a -io  -a "${breakpoint_file}" -b "${repeat_file}" > "${out_file}"
