

%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/Multiple_Thetacs/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/Multiple_Thetacs/';


%Specify the groups of firms that have different coefficients
F_group_cell = {1;2;3;4;5;6;7;8;9};


%num_F_groups = size(F_group_cell,1);
%l = [1; zeros(num_F_groups-1,1); mean_g];
xlim_graph = [-10;150];
lp_confidence_sets_script_multiple_thetacs;


%Do the confidence sets for theta_g only
lp_confidence_sets_scrips_multiple_thetacs_for_thetag;