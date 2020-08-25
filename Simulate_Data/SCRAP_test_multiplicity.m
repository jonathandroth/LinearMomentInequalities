l1 =[ones(num_F_groups_parameters,1) / num_F_groups_parameters; mean_g] ;
l2 =[zeros( num_F_groups_parameters ,1) ; 1];
for(superSim = 1:20)

superSim
T = 1000 + burnout;

if(onLaptop == 0)
    numSimulationsIDSet = 5000;
else
    numSimulationsIDSet = 50;
end

y_bar_cell = cell(numSimulationsIDSet,1);
X_bar_cell = cell(numSimulationsIDSet,1);


baseSeed = randi([1 10^8]);
parfor(sim = 1:numSimulationsIDSet)


[J_t_array, J_tminus1_array, Pi_array, F_array, G_array, Pi_star_array, Eta_jt_shocks_array, Eta_t_vec] = simulate_data(sim + baseSeed, F , J ,T, burnout, sigma_nu, sigma_eps, sigma_w, sigma_zetaj, sigma_zetajft, rho, lambda,theta_c,theta_g, g_vec, mu_f);

[~, ~, y_bar_s, X_bar_s] = lp_create_moments_for_identified_set_fn(F_group_cell_moments, F_group_cell_parameters, F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array, use_basic_moments, lambda, combine_theta_g_moments);

y_bar_cell{sim,1} = y_bar_s;
X_bar_cell{sim,1} = X_bar_s;

    
end

y_bar = cellReduce( y_bar_cell, @(x,d) mean(x,d) );
X_bar = cellReduce( X_bar_cell, @(x,d) mean(x,d) );


N = (T - burnout) * numSimulationsIDSet;
cutoff = log(N)/sqrt(N);

[ci1, slack_l1, slack_u1] = cs_linear_delta_lp_fn(y_bar,X_bar,l1,0);
[ci2, slack_l2, slack_u2] = cs_linear_delta_lp_fn(y_bar,X_bar,l2,0);

tol = 10^-4;

whichBinding1 = abs(slack_u1) < tol;
whichBinding2 = abs(slack_u2) < tol;
if(superSim == 1)
  bindingMat1 = whichBinding1;
  bindingMat2 = whichBinding2;
else
    bindingMat1 = [bindingMat1, whichBinding1];
    bindingMat2 = [bindingMat2, whichBinding2];
end

end

[mean(bindingMat1, 2),...
mean(bindingMat2, 2)]