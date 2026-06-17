#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=12
#SBATCH --time=02:30:00
#SBATCH --mem=45GB
#SBATCH --export=ALL
#SBATCH --job-name=ss

module load compiler/intel/
module load numlib/mkl

source ~/anaconda3/bin/activate

ncores_x_grid=1
ncores_liouvillian=12

export OMP_NUM_THREADS=$ncores_liouvillian
export MKL_NUM_THREADS=$ncores_liouvillian

ulimit -s unlimited

python3 friction_heom_main.py ss
python3 friction_heom_main.py markovian
