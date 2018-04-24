load( char(strcat( data_input_dir, dirname, 'all_simulation_params') ) ) %load the parameters needed to simulate the data


%Draw 5000 1000-period chains to compute identified set
T = 1000 + burnout;

if(onLaptop == 0)
    numSimulationsIDSet = 5000;
else
    numSimulationsIDSet = 2; %If on laptop, do smaller simulation
end

y_bar_cell = cell(numSimulationsIDSet,1);
X_bar_cell = cell(numSimulationsIDSet,length(lambda_vec));

parfor(sim = 1:numSimulationsIDSet )
%for(sim = 1:numSimulationsIDSet)

[J_t_array, J_tminus1_array, Pi_array, F_array, G_array, Pi_star_array, Eta_jt_shocks_array, Eta_t_vec] = simulate_data(sim, F , J ,T, burnout, sigma_nu, sigma_eps, sigma_w, sigma_zetaj, sigma_zetajft, rho, lambda_true,theta_c,theta_g, g_vec, mu_f);
    
    X_bar_cell_row = cell(1,length(lambda_vec) );
    
    for lambda_index = 1:length(lambda_vec)
        lambda = lambda_vec(lambda_index);
        
        [~, ~, y_bar_s, X_bar_s] = lp_create_moments_for_identified_set_fn(F_group_cell_moments, F_group_cell_parameters, F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array, use_basic_moments, lambda, combine_theta_g_moments);
        
        X_bar_cell_row{1,lambda_index} = X_bar_s;

        
    end
    
        y_bar_cell{sim,1} = y_bar_s;
        X_bar_cell(sim,:) = X_bar_cell_row;
end

y_bar = cellReduce( y_bar_cell, @(x,d) mean(x,d) );
X_bar_collapsed_cell = cell(1, length(lambda_vec));

for lambda_index = 1:length(lambda_vec)
    X_bar_collapsed_cell{1,lambda_index} = cellReduce( X_bar_cell(:,lambda_index), @(x,d) mean(x,d) );
end

N = (T - burnout) * numSimulationsIDSet;

for lambda_index = 1:length(lambda_vec)
    X_bar_lambda = X_bar_collapsed_cell{lambda_index};
    lambda_identified_set(lambda_index) = ( test_delta_lp_fn(y_bar, X_bar_lambda) <= (log(N)/sqrt(N)) );
    lambda_identified_set_zerocutoff(lambda_index) = ( test_delta_lp_fn(y_bar, X_bar_lambda) <= 0 );
end

save( strcat(ds_dir, 'lambda_identified_set'), ...
    'lambda_identified_set', 'lambda_identified_set_zerocutoff', 'lambda_vec');
