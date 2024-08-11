#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=2:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

module load python/3.10.9_anaconda2023.03_libmamba

source ~/software.install/py_venv/mlst/bin/activate

species=calbicans
ref_genome=SC5314_A21

while IFS= read -r line; do
    file="$line"
    python pubmlst_rest.py "$file" "$species" "$ref_genome"
    sleep 5
done < fasta_subsets.txt
