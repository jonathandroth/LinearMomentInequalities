#!/bin/bash 
#SBATCH -o ../../Logs/basic_inequalities3.out 
#SBATCH -e ../../Logs/basic_inequalities3.err 
#SBATCH -p serial_requeue
#SBATCH -c 24 
#SBATCH -N 1 
#SBATCH --mem 20000 
#SBATCH -t 2-01:00

source new-modules.sh
module load matlab
srun -n 1 -c 24 matlab-default -nodisplay -nosplash <run_basic_inequalities_main3.m
