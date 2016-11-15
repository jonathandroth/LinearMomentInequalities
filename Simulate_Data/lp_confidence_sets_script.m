basic_inequalities_set_parameters   

diagonal = 0;

working_dir = '/Users/jonathanroth/Google Drive/Research Projects/Moment_Inequalities_Ariel/Code/Simulate_Data';
%working_dir = '/n/home12/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data';
%parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

cd( working_dir);
%Specify where the input data is (can be relative to the working_dir)
data_input_dir = '../../Output/Simulated_Data/';

%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/';

dirname = 'Calibrated_SigmaZeta/';

%%
numdatasets = 1;
tic;


    
    
    long_ds_object = load( char(strcat( data_input_dir, dirname, 'ds_long.mat') )) ;
    length_long_chain = size( long_ds_object.J_t_array,3);
    
    
    
    
    %Create cells to store y_T, X_T, c_lf, and Sigma for each of the
    %datasets
    y_T_cell = cell(numdatasets,1);
    X_T_cell = cell(numdatasets,1);
    c_lf_cell = cell(numdatasets,1);
    Sigma_conditional_cell = cell(numdatasets,1);
    
    
    
    %Create cells in which to load the simulated dataset items
    F_array_cell = cell(numdatasets,1);
    G_array_cell = cell(numdatasets,1);
    Eta_jt_shocks_array_cell = cell(numdatasets,1);
    Eta_t_vec_cell = cell(numdatasets,1);
    J_t_array_cell = cell(numdatasets,1);
    J_tminus1_array_cell = cell(numdatasets,1);
    Pi_array_cell = cell(numdatasets,1);
    
    %The following loop takes subsets of the long chain and stores in a
    %cell. Each individual dataset will then be passed to a worker in the
    %parfor loop below
    for ds = 1:numdatasets

    rng(ds);
    rand_index = randsample( length_long_chain , nummarkets);
    
    F_array_cell{ds} = long_ds_object.F_array(:,:,rand_index);
    G_array_cell{ds} = long_ds_object.G_array(:,:,rand_index);
    Eta_jt_shocks_array_cell{ds} = long_ds_object.Eta_jt_shocks_array(:,:,rand_index);
    Eta_t_vec_cell{ds} = long_ds_object.Eta_t_vec(rand_index);
    %Pi_star_array = long_ds_object.Pi_star_array(:,:,rand_index);
    J_t_array_cell{ds} = long_ds_object.J_t_array(:,:,rand_index);
    J_tminus1_array_cell{ds} = long_ds_object.J_tminus1_array(:,:,rand_index);
    Pi_array_cell{ds} = long_ds_object.Pi_array(:,:,rand_index);
    end
    clear long_ds_object    
    
    
    parfor ds = 1:numdatasets
        
    F_array = F_array_cell{ds};
    G_array = G_array_cell{ds};
    Eta_jt_shocks_array = Eta_jt_shocks_array_cell{ds};
    Eta_t_vec = Eta_t_vec_cell{ds};
    %Pi_star_array = Pi_star_array_cell{ds};
    J_t_array = J_t_array_cell{ds};
    J_tminus1_array = J_tminus1_array_cell{ds};
    Pi_array = Pi_array_cell{ds};
   
    ds
    
    [moment_fn_allparams,moment_fn_interacted_allparams,Y, A_g, A_c,Y_basic, A_g_basic, A_c_basic] = ...
        generate_moment_fn( F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array); 
    moment_fn = @(theta_c, theta_g) -moment_fn_allparams(theta_c, theta_g, lambda_true);

   
        %A_g and A_c come out as functions of labmda
        %We replace these with their value at the true lambda
        A_g = A_g(lambda_true);
        A_c = A_c(lambda_true); 
        A_g_basic = A_g_basic(lambda_true);
        A_c_basic = A_c_basic(lambda_true); 
        
        
        %Merge theta_g and theta_c coefficients
        A = [A_c, A_g];
        A_basic = [A_c_basic, A_g_basic];
        
        %Remove any all zero columns
        A= A(:,any(A));
        A_basic = A_basic(:,any(A_basic));
        
        %Calculate the variance of Y conditional on A using the Abadie et
        %al matched pairs method
        Sigma_conditional = conditional_variance_fn(Y, A, diagonal);
        Sigma_conditional_basic  =conditional_variance_fn(Y_basic, A_basic, diagonal);
         
    
  %y_T and x_T are constructed so that population moments = y_T - X_T * delta
    % WE construct these so that they are less than 0 in expectation (the
    % moment fns are greater than 0 in expectation)
  T = size(A_g,1);
  X_T = [ sum( A_c ,1)', sum(A_g,1)'] / sqrt( T ); 
  y_T = -sum(Y,1)' / sqrt(T);
  %max( abs( (y_T + X_T * [1 ;2]) - mean(  moment_fn_interacted_allparams(1,2, lambda_true) )' ))
    
  X_T_cell(ds,1) = {X_T};
  y_T_cell(ds,1) = {y_T};
  Sigma_conditional_cell(ds,1) = {Sigma_conditional};
    %     
%     ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, 'Interacted_Moments/grid', num2str(ds));
%     ds_name = ds_name{:};
%     save( ds_name, 'grid_lf', 'grid_rsw', 'grid_conditional', 'grid_hybrid');
    end

    %Save the cells for the interacted moments

    ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/values_for_lp');
    save( ds_name, 'y_T_cell', 'X_T_cell', 'Sigma_conditional_cell', 'c_lf_cell');
    


toc;


%% 
mean_g = mean(G_array(1,:,1))
l = [1; mean_g];
delta_true = [theta_c_true; theta_g_true];
l'*delta_true

for ds = 1:numdatasets
    
   X_T = X_T_cell{ds,1};
   y_T = y_T_cell{ds,1};
   Sigma = Sigma_conditional_cell{ds,1}; 
   
   c_alpha = c_lf(Sigma, alpha, Z_draws_interacted);
   
   test_delta_lp_fn(y_T, X_T);
   
 
   cs_linear_delta_lp_fn(y_T,X_T,l,c_alpha);
   
   c_lp_alpha = c_lf_lp(X_T,Z_draws_interacted(:,1:200),Sigma,alpha);
   cs_linear_delta_lp_fn(y_T,X_T,l,c_lp_alpha);
end
