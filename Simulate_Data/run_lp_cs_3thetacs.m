

%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/3Thetacs/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/3Thetacs/';


%Specify the groups of firms that have different coefficients
F_group_cell = {[1;2;3];[4;5;6];[7;8;9]};


%num_F_groups = size(F_group_cell,1);
%l = [1; zeros(num_F_groups-1,1); mean_g];

lp_confidence_sets_script_multiple_thetacs;