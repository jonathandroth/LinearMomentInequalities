#!/bin/bash 
#SBATCH -o ../../Logs/lp_cs_multiple_thetas_combined_thetag.out
#SBATCH -e ../../Logs/lp_cs_multiple_thetas_combined_thetag.err
#SBATCH -p shared
#SBATCH -c 24 
#SBATCH -N 1 
#SBATCH --mem 128000 
#SBATCH -t 2-01:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jonathanroth@g.harvard.edu

source new-modules.sh
module load matlab
srun -n 1 -c 24 matlab -nodisplay -nosplash <run_lp_cs_multiple_thetas_combined_thetag.m
