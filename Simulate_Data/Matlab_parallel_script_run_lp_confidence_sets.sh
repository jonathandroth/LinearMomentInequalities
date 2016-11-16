#!/bin/bash 
#SBATCH -o ../../Logs/lp_confidence_sets.out
#SBATCH -e ../../Logs/lp_confidence_sets.err
#SBATCH -p serial_requeue
#SBATCH -c 24 
#SBATCH -N 1 
#SBATCH --mem 20000 
#SBATCH -t 2-01:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jonathanroth@g.harvard.edu

source new-modules.sh
module load matlab
srun -n 1 -c 24 matlab-default -nodisplay -nosplash <lp_confidence_sets_script.m
