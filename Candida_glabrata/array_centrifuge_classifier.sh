#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -p msismall,msilarge
#SBATCH --array=

centrifuge_home=/home/selmecki/shared/software/centrifuge
db=/home/selmecki/shared/centrifuge_fungi_refseq/ncbi_nt/nt
line=${SLURM_ARRAY_TASK_ID}
sample_file=

#Read sample file line corresponding to array task ID and get variables
strain=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)

srun "${centrifuge_home}"/centrifuge -p 10 -x "${db}" \
-1 ../trimmed_fastq/"${strain}"_trimmed_1P.fq \
-2 ../trimmed_fastq/"${strain}"_trimmed_2P.fq \
-S "${strain}"_centrifuge.txt

srun "${centrifuge_home}"/centrifuge-kreport -x "${db}" "${strain}"_centrifuge.txt > "${strain}"_kreport.txt
