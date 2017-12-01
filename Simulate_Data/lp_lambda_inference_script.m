%% Specify parameters
lp_set_parameters
%% Import data and calculate moments and covariance matrices

%lambda_vec = (0.01:10:20.01)';
%lambda_vec = (0.00:0.01:0.01)';
%lambda_vec = (0.01:10:10.01)';
%lambda_vec = (0.01:.2:5.01)';

conditional_test = NaN(numdatasets, size(lambda_vec,1));
hybrid_test = NaN(size(conditional_test));
lf_test_original = NaN(size(conditional_test));
lf_test_modified = NaN(size(conditional_test));

lambda_identified_set = NaN( size(lambda_vec) );
lambda_identified_set_zerocutoff = NaN( size(lambda_vec) );

i = 1;


for lambda = lambda_vec'
%lambda_true = .386

lambda
lp_create_moments_and_covariances

nummoments = size( y_T_cell{ds,1} , 1);
rng(0);
Z_draws_interacted = randn(nummoments, 10000);

parfor ds = 1:numdatasets
X_T = X_T_cell{ds};
y_T = y_T_cell{ds};
Sigma = Sigma_conditional_cell{ds};


conditional_test(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);
hybrid_test(ds,i) = lp_hybrid_test_fn( y_T, X_T, Sigma, alpha, alpha/10);

%Do non-conditional tests
c_alpha = c_lf(Sigma, alpha, Z_draws_interacted);
c_lp_alpha = c_lf_lp(X_T,Z_draws_interacted(:,1:numsims_lp),Sigma,alpha);
eta = test_delta_lp_fn( y_T, X_T);

lf_test_original(ds,i) = eta > c_alpha;
lf_test_modified(ds,i) = eta > c_lp_alpha;

end

i = i+1;
end

%% Save the results

    ds_dir = strcat( data_output_dir, 'Interacted_Moments/');
    mkdir(ds_dir);
    save( strcat(ds_dir, 'lambda_results'), 'conditional_test', 'lf_test_original', 'lf_test_modified', 'hybrid_test', 'lambda_vec');

    

%% Compute identified set

tic;
load( char(strcat( data_input_dir, dirname, 'all_simulation_params') ) ) %load the parameters needed to simulate the data


%Draw 5000 1000-period chains to compute identified set
T = 1000 + burnout;
%numSimulations = 5000; %CHANGE ME BACK
numSimulations = 50;

y_bar_cell = cell(numSimulations,1);
X_bar_cell = cell(numSimulations,length(lambda_vec));

parfor(sim = 1:numSimulations)
%for(sim = 1:numSimulations)

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

N = (T - burnout) * numSimulations;

for lambda_index = 1:length(lambda_vec)
    X_bar_lambda = X_bar_collapsed_cell{lambda_index};
    lambda_identified_set(lambda_index) = ( test_delta_lp_fn(y_bar, X_bar_lambda) <= (log(N)/sqrt(N)) );
    lambda_identified_set_zerocutoff(lambda_index) = ( test_delta_lp_fn(y_bar, X_bar_lambda) <= 0 );
end

save( strcat(ds_dir, 'lambda_identified_set'), ...
    'lambda_identified_set', 'lambda_identified_set_zerocutoff', 'lambda_vec');
    
    
display('Finished identified set calc');

toc;


 %% Plot

 ds_dir = strcat( data_output_dir, 'Interacted_Moments/');
 load(  strcat(ds_dir, 'lambda_identified_set') );
 load(  strcat(ds_dir, 'lambda_results') );
 
  
 plot(lambda_vec, [mean(conditional_test);...
                  mean(hybrid_test);...
                  mean(lf_test_original);...
                  mean(lf_test_modified)] ) 


identified_set_max = max( lambda_vec( lambda_identified_set == 1) );
identified_set_min = min( lambda_vec( lambda_identified_set == 1) );


if(~isempty(identified_set_max) )
    line( [identified_set_max; identified_set_max], [0;1], 'LineStyle', '--', 'Color',  'r');
else
    warning('Didnt find any lambdas in identfied set');
end

if(~isempty(identified_set_min) && identified_set_min ~= min(lambda_vec)  )
    line( [identified_set_min; identified_set_min], [0;1], 'LineStyle', '--', 'Color',  'r');
end

legend( 'Conditional', 'Hybrid', 'LF','LFN', 'Identified Set Bound', 'Location','eastoutside' );
ylabel('Rejection Probability');
xlabel('Lambda');
title('Rejection Probabilities for Lambda');
saveas( gcf, strcat(figures_output_dir,filename_graph ), 'epsc');

display( strcat(figures_output_dir,filename_graph ) )
