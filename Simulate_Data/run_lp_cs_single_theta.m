

%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/Single_Thetac/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/Single_Thetac/';


%Specify the groups of firms that have different coefficients
F_group_cell = {[1;2;3;4;5;6;7;8;9]};


%num_F_groups = size(F_group_cell,1);
%l = [1; zeros(num_F_groups-1,1); mean_g];

%Run the main script (this does confidence sets for the mean weight)
lp_confidence_sets_script_multiple_thetacs;


%Do the confidence sets for theta_g only
l = [zeros( size(F_group_cell,1) ,1) ; 1];
xlim_graph = [-150;100];
filename_graph =  'Theta_g Rejection Probabilities';
xlabel_graph = 'Theta_g';

lp_compute_confidence_sets;
lp_graphs;