#!/bin/bash 
#SBATCH -o ../../Logs/run_lp_cs_3thetacs_interactedmoments.out
#SBATCH -e ../../Logs/run_lp_cs_3thetacs_interactedmoments.err
#SBATCH -p shared
#SBATCH -c 24 
#SBATCH -N 1 
#SBATCH --mem 80000 
#SBATCH -t 2-01:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jonathanroth@g.harvard.edu

source new-modules.sh
module load matlab
srun -n 1 -c 24 matlab-default -nodisplay -nosplash <run_lp_cs_3thetacs_interactedmoments.m
