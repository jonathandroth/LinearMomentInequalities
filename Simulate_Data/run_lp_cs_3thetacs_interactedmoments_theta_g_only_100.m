%offset for seed 
dsoffset = 100;

% Toggle doing only AS/KMS
as_kms_only = 1;

%specify whether to combine the moments for theta_g, or to have separate
%ones for each firm
combine_theta_g_moments = 1;

%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/3Thetacs/Interacted_Moments/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/3Thetacs/Interacted_Moments/';


%Specify the groups of firms that have different coefficients
F_group_cell_moments = {[1;2;3];[4;5;6];[7;8;9]};


%num_F_groups_moments = size(F_group_cell_moments,1);
%l = [1; zeros(num_F_groups_moments-1,1); mean_g];

if( exist('xlim_graph_meanweight'))
    xlim_graph = xlim_graph_meanweight;
else
    xlim_graph = [-150;175];
end

if( exist('xsplit_graph_meanweight'))
   xsplit_graph = xsplit_graph_meanweight;
end

filename_graph = 'Mean_Weight_Rejection_Probabilities';
%lp_confidence_sets_script_multiple_thetacs;


%Do the confidence sets for theta_g only
lp_confidence_sets_script_multiple_thetacs_for_thetag_only;
