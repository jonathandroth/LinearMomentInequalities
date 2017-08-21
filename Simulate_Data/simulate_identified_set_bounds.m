%1thetac_interacted params

%specify whether to combine the moments for theta_g, or to have separate
%ones for each firm
combine_theta_g_moments = 1;


%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/1Thetac/Interacted_Moments/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/1Thetac/Interacted_Moments/';


%Specify the groups of firms that have different coefficients
F_group_cell_moments = {[1;2;3;4;5;6;7;8;9]};


%num_F_groups_moments = size(F_group_cell_moments,1);
%l = [1; zeros(num_F_groups_moments-1,1); mean_g];

%Run the main script (this does confidence sets for the mean weight)
if( exist('xlim_graph_meanweight'))
    xlim_graph = xlim_graph_meanweight;
else
    xlim_graph = [-150;175];
end

if( exist('xsplit_graph_meanweight'))
   xsplit_graph = xsplit_graph_meanweight;
end

filename_graph = 'Mean_Weight_Rejection_Probabilities';

lp_set_parameters
%%
load('../../Output/Simulated_Data/all_simulation_params')


numSimulations = 2;

boundsMat_meanweight = NaN(numSimulations,2);
ub_slackMatrix_meanweight = [];
lb_slackMatrix_meanweight = [];
boundsMat_thetag = NaN(numSimulations,2);
ub_slackMatrix_thetag = [];
lb_slackMatrix_thetag = [];


for(sim = 1:numSimulations)
[J_t_array, J_tminus1_array, Pi_array, F_array, G_array, Pi_star_array, Eta_jt_shocks_array, Eta_t_vec] = simulate_data(1, F , J ,T, burnout, sigma_nu, sigma_eps, sigma_w, sigma_zetaj, sigma_zetajft, rho, lambda,theta_c,theta_g, g_vec, mu_f);

lp_create_moments_for_identified_set;

%Compute ID set for meanweight
mean_g = mean(G_array(1,:,1));
l_meanweight = [ones(num_F_groups_parameters,1) / num_F_groups_parameters; mean_g];
[identified_set_bounds_for_meanweight, slack_lb_meanweight, slack_ub_meanweight] = cs_linear_delta_lp_fn(y_T,X_T,l_meanweight,0);


boundsMat_meanweight(sim,:) = identified_set_bounds_for_meanweight;
ub_slackMatrix_meanweight = [ ub_slackMatrix_meanweight , slack_ub_meanweight];
lb_slackMatrix_meanweight = [ lb_slackMatrix_meanweight , slack_lb_meanweight];


l_thetag = [ zeros( num_F_groups_parameters ,1) ; 1 ];
[identified_set_bounds_for_thetag, slack_lb_thetag, slack_ub_thetag] = cs_linear_delta_lp_fn(y_T,X_T,l_thetag,0);


boundsMat_thetag(sim,:) = identified_set_bounds_for_thetag;
ub_slackMatrix_thetag = [ ub_slackMatrix_thetag , slack_ub_thetag];
lb_slackMatrix_thetag = [ lb_slackMatrix_thetag , slack_lb_thetag];

end


