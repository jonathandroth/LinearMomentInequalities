#!/bin/bash 
#SBATCH -o ../../Logs/differences_grids.out 
#SBATCH -e ../../Logs/differenes_grids.err 
#SBATCH -p shared
#SBATCH -c 1 
#SBATCH -N 1 
#SBATCH --mem 20000 
#SBATCH -t 1-01:00

source new-modules.sh
module load matlab
srun -n 1 -c 1  matlab-default -nodisplay -nosplash < rejection_grids_differences.m
