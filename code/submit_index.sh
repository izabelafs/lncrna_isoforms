#!/bin/bash -l
#SBATCH -J index
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=izabela.silva@uni.lu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=30GB
#SBATCH --time=1-00:00:00
#SBATCH -p batch

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"
# Your more useful application can be started below!

conda activate salmon
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
conda deactivate

