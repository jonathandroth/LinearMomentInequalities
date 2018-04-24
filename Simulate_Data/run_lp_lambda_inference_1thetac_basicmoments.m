%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/1Thetac/Basic_Moments/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/1Thetac/Basic_Moments/';

filename_graph = 'lambda_rejection_probabilities';

%Specify the groups of firms that have different coefficients
F_group_cell_moments = {[1;2;3;4;5;6;7;8;9]};


%num_F_groups_moments = size(F_group_cell_moments,1);
%l = [1; zeros(num_F_groups_moments-1,1); mean_g];

use_basic_moments = 1;
combine_theta_g_moments = 1;

lambda_vec = linspace(15,140,100)';
lp_lambda_inference_script