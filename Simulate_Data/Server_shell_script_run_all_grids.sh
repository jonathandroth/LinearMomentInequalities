#!/bin/bash
# This script initiates on the Odyssey server the various permutations of the grid search over two parameters

#Run using Full conditional covariance matrix
sbatch Matlab_parallel_script_run_basic_inequalities.sh

#Run using unconditional covariance matrix
sbatch Matlab_parallel_script_run_basic_inequalities2.sh

#Run using a conditional covariance matrix where the Malahanobis distance only depends on the diagonalized covariance matrix of the X's
sbatch Matlab_parallel_script_run_basic_inequalities3.sh

#Run using the "Oracle" conditional covariance matrix
sbatch Matlab_parallel_script_run_basic_inequalities4.sh

