#!/bin/bash
# This script initiates on the Odyssey server the various permutations of the linear programming methods

# Inference holding lambda constant

#1 thetac
sbatch Matlab_parallel_script_run_lp_cs_1thetac_interactedmoments.sh
sbatch Matlab_parallel_script_run_lp_cs_1thetac_basicmoments.sh

#3 thetacs
sbatch Matlab_parallel_script_run_lp_cs_3thetacs_interactedmoments.sh
sbatch Matlab_parallel_script_run_lp_cs_3thetacs_basicmoments.sh

#9 thetacs
sbatch Matlab_parallel_script_run_lp_cs_9thetacs_interactedmoments.sh
sbatch Matlab_parallel_script_run_lp_cs_9thetacs_basicmoments.sh


# Inference with lambda
sbatch Matlab_parallel_script_run_lp_lambda_inference_single_theta.sh
sbatch Matlab_parallel_script_run_lp_lambda_inference_3thetacs.sh
sbatch Matlab_parallel_script_run_lp_lambda_inference_9thetacs.sh
