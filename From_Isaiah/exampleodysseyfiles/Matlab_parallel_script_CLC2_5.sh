#!/bin/bash 
#SBATCH -o calc2_5.out 
#SBATCH -o calc2_5.err 
#SBATCH -p serial_requeue
#SBATCH -n 24 
#SBATCH -N 1 
#SBATCH --mem 20000 
#SBATCH -t 7-00:00

module load math/matlab-R2014b
matlab-default -nodisplay -nosplash <CLC_CS_Calc_parallel_5points_2.m