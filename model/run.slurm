#!/bin/bash
#SBATCH --job-name=training
#SBATCH --account=stf
#SBATCH --partition=gpu-2080ti
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=01:00:00
#SBATCH --mem=150GB
#SBATCH -o /gscratch/scrubbed/wgalvin/model/train.out
# You can get email notifications when your job starts and stops - useful for long running jobs
#SBATCH --mail-user=wgalvin@uw.edu
#SBATCH --mail-type=ALL

python protein_holography_pytorch/training/chi_training.py

exit 0