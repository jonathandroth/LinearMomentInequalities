
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
l_meanweight = [ones(num_F_groups_parameters,1) / num_F_groups_parameters; mean_g];
[identified_set_bounds_for_meanweight, slack_lb_meanweight, slack_ub_meanweight] = cs_linear_delta_lp_fn(y_T,X_T,l_meanweight,0)';


boundsMat_meanweight(sim,:) = identified_set_bounds_for_meanweight;
ub_slack_matrix_meanweight = [ ub_slack_matrix_meanweight , slack_ub_meanweight];
lb_slack_matrix_meanweight = [ lb_slack_matrix_meanweight , slack_lb_meanweight];


l_thetag = [ zeros( num_F_groups_parameters ,1) ; 1 ];
[identified_set_bounds_for_thetag, slack_lb_thetag, slack_ub_thetag] = cs_linear_delta_lp_fn(y_T,X_T,l_thetag,0)';


boundsMat_thetag(sim,:) = identified_set_bounds_for_thetag;
ub_slack_matrix_thetag = [ ub_slack_matrix_thetag , slack_ub_thetag];
lb_slack_matrix_thetag = [ lb_slack_matrix_thetag , slack_lb_thetag];

end


