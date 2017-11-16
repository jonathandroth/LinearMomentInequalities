#!/bin/bash 
#SBATCH -o ../../Logs/basic_inequalities2.out 
#SBATCH -e ../../Logs/basic_inequalities2.err 
#SBATCH -p shared
#SBATCH -c 24 
#SBATCH -N 1 
#SBATCH --mem 40000 
#SBATCH -t 2-01:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jonathanroth@g.harvard.edu

source new-modules.sh
module load matlab
srun -n 1 -c 24  matlab-default -nodisplay -nosplash <run_basic_inequalities_main2.m
