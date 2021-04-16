%This script computes the confidence sets and identified set for l * theta
%in the lp_confidence_sets script

%If l is not specfied in the environment, it is assumed to be [1/F,...,1/F,
%mean_g]

%This script used to be part of lp_confidence_sets_multiple_thetacs, but I
%broke it out so that we can run it multiple times with different l's,
%after loading the data

addpath(genpath(pwd))

G_array = G_array_cell{1};
mean_g = mean(G_array(1,:,1));

%If l not specified, then do the l that gives you the mean weights
if( exist('l') == 0)
l = [ones(num_F_groups_parameters,1) / num_F_groups_parameters; mean_g];
%l = [1;0];
end

delta_true = [repmat(theta_c_true,num_F_groups_parameters,1);  theta_g_true];
l'*delta_true


%% Compute identified set

tic;
load( char(strcat( data_input_dir, dirname, 'all_simulation_params') ) ) %load the parameters needed to simulate the data


%Draw 5000 1000-period chains to compute identified set
T = 1000 + burnout;

if(onLaptop == 0)
    numSimulationsIDSet = 5000;
else
    numSimulationsIDSet = 2;
end

y_bar_cell = cell(numSimulationsIDSet,1);
X_bar_cell = cell(numSimulationsIDSet,1);

parfor(sim = 1:numSimulationsIDSet)


[J_t_array, J_tminus1_array, Pi_array, F_array, G_array, Pi_star_array, Eta_jt_shocks_array, Eta_t_vec] = simulate_data(sim, F , J ,T, burnout, sigma_nu, sigma_eps, sigma_w, sigma_zetaj, sigma_zetajft, rho, lambda,theta_c,theta_g, g_vec, mu_f);

[~, ~, y_bar_s, X_bar_s] = lp_create_moments_for_identified_set_fn(F_group_cell_moments, F_group_cell_parameters, F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array, use_basic_moments, lambda, combine_theta_g_moments);

y_bar_cell{sim,1} = y_bar_s;
X_bar_cell{sim,1} = X_bar_s;

    
end

y_bar = cellReduce( y_bar_cell, @(x,d) mean(x,d) );
X_bar = cellReduce( X_bar_cell, @(x,d) mean(x,d) );


N = (T - burnout) * numSimulationsIDSet;
cutoff = log(N)/sqrt(N);

%Compute ID set bounds with log(N)/sqrt(N) cutoff
identified_set_bounds = cs_linear_delta_lp_fn(y_bar,X_bar,l,cutoff)';
ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/identified_set_bounds');
    mkdir(ds_name);
    save( ds_name, 'identified_set_bounds');

%Compute ID set bounds with 0 cutoff
identified_set_bounds_zerocutoff = cs_linear_delta_lp_fn(y_bar,X_bar,l,0)';
ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/identified_set_bounds_zerocutoff');
    mkdir(ds_name);
    save( ds_name, 'identified_set_bounds_zerocutoff');
    
    
display('Finished identified set calc');

toc;
%% Confidence sets for linear combination of theta's

%load( strcat( data_output_dir, dirname, 'Interacted_Moments/values_for_lp') )
%load( strcat( data_output_dir, dirname, 'Interacted_Moments/grid_cell') )

%For the LF methods, we save bounds of confidence sets, since we compute
%these exactly. These are intiated here
confidence_sets_using_c_alpha = NaN(numdatasets,2);
confidence_sets_using_c_lp_alpha = NaN(numdatasets,2);

%We do the same thing for the projection methods using the KMS package
confidence_sets_using_as = NaN(numdatasets,2);
confidence_sets_using_kms = NaN(numdatasets,2);


%For the conditional/hybrid methods, we get rejection probabilities over a
%grid of values for l*theta. So we initialize those grids here
num_beta0_gridpoints = 1001;
beta0_grid = linspace( xlim_graph(1), xlim_graph(2), num_beta0_gridpoints);

%Add ID set endpoints using both methods to the grid
beta0_grid = sort([beta0_grid, identified_set_bounds, identified_set_bounds_zerocutoff]); 
num_beta0_gridpoints = length(beta0_grid); %update length for the id set endpoints

rejection_grid_conditional = NaN(numdatasets, num_beta0_gridpoints);
rejection_grid_hybrid = NaN(numdatasets, num_beta0_gridpoints);

rejection_grid_rcc = NaN(numdatasets, num_beta0_gridpoints);
rejection_grid_cc = NaN(numdatasets, num_beta0_gridpoints);


nummoments = size( y_T_cell{ds,1} , 1);
Z_draws_interacted = randn(nummoments, 10000);

%%CHANGE THIS BACK WHEN DONE DEBUGGING
%for ds = 1:numdatasets
parfor ds = 1:numdatasets

ds
   X_T = X_T_cell{ds,1};
   y_T = y_T_cell{ds,1};
   Sigma = Sigma_conditional_cell{ds,1}; 
   
   %Rescale moments by the variance (this is done earlier now)
%    D_sigma_minushalf = diag( sqrt( diag(Sigma) ).^(-1) );
%    y_T = D_sigma_minushalf * y_T;
%    X_T = D_sigma_minushalf * X_T;
%    
   
   %Do the LF (aka LFP) test and store the resulting confidence set
   c_alpha = c_lf(Sigma, alpha, Z_draws_interacted); 
   confidence_sets_using_c_alpha(ds,:) = cs_linear_delta_lp_fn(y_T,X_T,l,c_alpha)';
   
   %Do the modified least favorable test (aka LF) and store the resulting confidence
   %set
   c_lp_alpha = c_lf_lp(X_T,Z_draws_interacted(:,1:numsims_lp),Sigma,alpha);
   confidence_sets_using_c_lp_alpha(ds,:) = cs_linear_delta_lp_fn(y_T,X_T,l,c_lp_alpha)';
      
 

   %For the KMS/AMS and cond'l/hybrid, 
   %we do a change of basis such that M*theta yields l*theta in the first
   %row (I construct an orthonormal basis for the rest of the basis
   %vectors, although this is not strictly necessary)
   M = [l' ; null( repmat(l',size(l)) )'];
   X_T = X_T * M^(-1);
   
   
   %Do AS and KMS for specs with <9 parameters
   if size(X_T,2) < 9
   try
       [as_ci,as_output] = projected_AS_or_KMS(y_T, X_T, 500, Sigma,[1;zeros(size(X_T,2)-1,1)], NaN, 'AS');
       
       %If have a convergence error for either bound, set to NaN
       if (as_output.flagL_EAM ~= 1) || (as_output.flagU_EAM ~= 1)
           as_ci = [NaN, NaN];
       end
       confidence_sets_using_as(ds,:) = as_ci;
   catch
   end
   
   try
   
       [kms_ci,kms_output] = projected_AS_or_KMS(y_T, X_T, 500, Sigma,[1;zeros(size(X_T,2)-1,1)], NaN, 'KMS');
       
       %If have a convergence error for either bound, set to NaN
       if (kms_output.flagL_EAM ~= 1) || (kms_output.flagU_EAM ~= 1)
           kms_ci = [NaN, NaN];
       end
       confidence_sets_using_kms(ds,:) = kms_ci;
   catch
   end
   
   end
   %%%Do the conditional and hybrid tests treating l as a non-linear
   %parameter
    
    conditional_rejection_vec = NaN(num_beta0_gridpoints,1);
    hybrid_rejection_vec = NaN(num_beta0_gridpoints,1);
    rcc_rejection_vec = NaN(num_beta0_gridpoints,1);
    cc_rejection_vec = NaN(num_beta0_gridpoints,1);

    
    %Prior to looping through all of the value for l*theta, compute the
    %least-favorable critical values, since this doens't depend on the
    %value of l*theta. These will be inputted to the hybrid method
    X_T_tilde = X_T(:,2:end);
    [~, lf_simulated_draws] = c_lf_lp(X_T_tilde,Z_draws_interacted(:,1:numsims_lp),Sigma,alpha);
   
    %Loop over grid for beta_0 = l*theta values
    %For each beta_0, create y_Tilde = y_T - X_t[:,1] * b_0, x_Tilde =
    %X_T[:,2:]; and run the conditional/hybrid approaches for inference on
    %beta_0
    count = 1;
    for beta0 = beta0_grid
    
        y_T_tilde = y_T - X_T(:,1) * beta0;
        
        
        conditional_rejection_vec(count,1) = lp_conditional_test_fn( y_T_tilde, X_T_tilde, Sigma, alpha);
        hybrid_rejection_vec(count,1) = lp_hybrid_test_fn( y_T_tilde, X_T_tilde, Sigma, alpha, alpha/10, lf_simulated_draws);
        
        %[T_CC, c_RCC, c_CC] = rcc_test_fn(sqrt(nummarkets)^(-1) * y_T_tilde, eye(size(y_T_tilde,1)), X_T_tilde, zeros(size(y_T_tilde,1),1) , Sigma, nummarkets, 0, alpha);

        %If fewer than 10 params or fewer than 100 moments, do RCC test
            %Otherwise do CC test only and use ub of RCC test based on best
            %possible refinement
        if size(X_T,2) < 10 || size(X_T,1) < 100
            try
            [T_RCC,cv_RCC,cv_CC,~,~, dof_n] = func_subRCC(X_T_tilde, -y_T_tilde, Sigma, alpha);
            cc_rejection_vec(count,1) = (dof_n > 0) * (T_RCC > cv_CC);
            rcc_rejection_vec(count,1) = (dof_n > 0) * (T_RCC > cv_RCC);
            catch
            [T_CC,cv_CC, dof_n] = func_subCC(X_T_tilde, -y_T_tilde, Sigma, alpha);
            cc_rejection_vec(count,1) = (dof_n > 0) * (T_CC > cv_CC);
            rcc_rejection_vec(count,1) = ((dof_n > 0) * (T_CC > cv_CC)) | ( (dof_n ==1) & (T_CC > chi2inv(1-2*alpha,dof_n)) ); %set RCC to 1 if dof_n == 1 and T_CC is above the 1-2*alpha cv    
            end
        else
            [T_CC,cv_CC, dof_n] = func_subCC(X_T_tilde, -y_T_tilde, Sigma, alpha);
            cc_rejection_vec(count,1) = (dof_n > 0) * (T_CC > cv_CC);
            rcc_rejection_vec(count,1) = ((dof_n > 0) * (T_CC > cv_CC)) | ( (dof_n ==1) & (T_CC > chi2inv(1-2*alpha,dof_n)) ); %set RCC to 1 if dof_n == 1 and T_CC is above the 1-2*alpha cv
        end
        
        count = count +1;
    end
    
    rejection_grid_conditional(ds,:) = conditional_rejection_vec;
    rejection_grid_hybrid(ds,:) = hybrid_rejection_vec;
    
    rejection_grid_rcc(ds,:) = rcc_rejection_vec;
    rejection_grid_cc(ds,:) = cc_rejection_vec;
    
end

%Create the confidence sets using the grid search
%grid_min_max
    
    %Average the rejection probabilities for the hybrid and conditional
    %methods
    full_rejection_grid_conditional = rejection_grid_conditional;
    full_rejection_grid_hybrid = rejection_grid_hybrid;
    
    rejection_grid_conditional = mean(rejection_grid_conditional,1);
    rejection_grid_hybrid = mean(rejection_grid_hybrid,1);
    
    
    full_rejection_grid_rcc = rejection_grid_rcc;
    full_rejection_grid_cc = rejection_grid_cc;
    
    rejection_grid_rcc = mean(rejection_grid_rcc,1);
    rejection_grid_cc = mean(rejection_grid_cc,1);
    

    ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/confidence_sets_lp');
    mkdir(ds_name);
    save( ds_name, 'confidence_sets_using_c_alpha', 'confidence_sets_using_c_lp_alpha',...
                   'rejection_grid_hybrid', 'rejection_grid_conditional','beta0_grid' ,...
                   'full_rejection_grid_conditional','full_rejection_grid_hybrid',...
                   'rejection_grid_rcc', 'rejection_grid_cc', ...
                   'full_rejection_grid_rcc', 'full_rejection_grid_cc',...
                   'confidence_sets_using_as', 'confidence_sets_using_kms');

    
    
    
