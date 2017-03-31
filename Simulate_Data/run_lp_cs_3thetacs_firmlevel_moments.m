

%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/3Thetacs/Firm_Level_Moments/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/3Thetacs/Firm_Level_Moments/';


%Specify the groups of firms that have different coefficients and moments
    %This says that there will moments separately for each firm, but only
    %one theta_c parameter
F_group_cell_moments = {1;2;3;4;5;6;7;8;9};
F_group_cell_parameters = {[1;2;3];[4;5;6];[7;8;9]};

%num_F_groups_moments = size(F_group_cell_moments,1);
%l = [1; zeros(num_F_groups_moments-1,1); mean_g];

%Run the main script (this does confidence sets for the mean weight)
xlim_graph = [-10;150];
lp_confidence_sets_script_multiple_thetacs;


%Do the confidence sets for theta_g only
lp_confidence_sets_script_multiple_thetacs_for_thetag;

