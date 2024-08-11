#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=2:00:00
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH -p amdsmall,amdlarge,amd512
#SBATCH --array=1-32

set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID}
sample_file=PRJNA361477_SRR_Acc_List.txt

#Load modules
module load sratoolkit

accession=$(awk -v val="$line" 'NR == val { print $0}' $sample_file)

fastq-dump --split-3 "${accession}"
