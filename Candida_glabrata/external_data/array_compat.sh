#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-234

set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID}
module load compatibility/mesabi-centos7

compat-exec ./wrapped_subsample.sh $line

