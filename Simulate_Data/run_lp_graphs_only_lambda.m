

graphs_only = 1
remote_to_server = 1


run_lp_lambda_inference_1thetac_basicmoments
run_lp_lambda_inference_3thetacs_basicmoments
run_lp_lambda_inference_9thetacs_basicmoments

run_lp_lambda_inference_1thetac_interactedmoments
run_lp_lambda_inference_3thetacs_interactedmoments
run_lp_lambda_inference_9thetacs_interactedmoments
% 
% make_power_table
use_zero_cutoff =0;
make_excess_length_table_lambda;
make_size_table_lambda;

use_zero_cutoff =1;
make_excess_length_table_lambda
make_size_table_lambda
