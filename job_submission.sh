#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=6
#SBATCH --time=30:00:00
#SBATCH --mem=50GB
#SBATCH --export=ALL
#SBATCH --job-name=friction_HEOM

module load compiler/intel/
module load numlib/mkl

source ~/anaconda3_older/bin/activate

ncores_q_grid=1
ncores_liouvillian=6

export OMP_NUM_THREADS=$ncores_q_grid
export MKL_NUM_THREADS=$ncores_liouvillian

ulimit -s 400000

python3 friction_heom_main.py

