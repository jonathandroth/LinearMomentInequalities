#!/bin/bash 
#SBATCH -o ../../Logs/run_lp_cs_3thetacs_interactedmoments_thetag_only_200.out
#SBATCH -e ../../Logs/run_lp_cs_3thetacs_interactedmoments_thetag_only_200.err
#SBATCH -p shared
#SBATCH -c 36 
#SBATCH -N 1 
#SBATCH --mem 80000 
#SBATCH -t 7-00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jonathanroth@g.harvard.edu

source new-modules.sh
module load matlab
srun -n 1 -c 36 matlab -nodisplay -nosplash <run_lp_cs_3thetacs_interactedmoments_theta_g_only_200.m
