%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/3Thetacs/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/3Thetacs/';

filename_graph = 'lambda_rejection_probabilities';

%Specify the groups of firms that have different coefficients
F_group_cell_moments = {[1;2;3];[4;5;6];[7;8;9]};

use_basic_moments = 0;

lp_lambda_inference_script
