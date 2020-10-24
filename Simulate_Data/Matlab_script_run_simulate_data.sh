#!/bin/bash 
#SBATCH -o ../../Logs/simulate_data.out 
#SBATCH -e ../../Logs/simulate_data.err 
#SBATCH -p shared
#SBATCH -c 1 
#SBATCH -N 1 
#SBATCH --mem 16000 
#SBATCH -t 2-01:00

source new-modules.sh
module load matlab
srun -n 1 -c 1 matlab -nodisplay -nosplash <run_simulate_data.m
