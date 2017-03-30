%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/Single_Thetac/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/Single_Thetac/';

filename_graph = 'lambda_rejection_probabilities';

%Specify the groups of firms that have different coefficients
F_group_cell_moments = {[1;2;3;4;5;6;7;8;9]};


%num_F_groups = size(F_group_cell_moments,1);
%l = [1; zeros(num_F_groups-1,1); mean_g];

use_basic_moments = 0;

lp_lambda_inference_script
