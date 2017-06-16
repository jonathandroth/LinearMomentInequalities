%This script computes the confidence sets and identified set for l * theta
%in the lp_confidence_sets script

%If l is not specfied in the environment, it is assumed to be [1/F,...,1/F,
%mean_g]

%This script used to be part of lp_confidence_sets_multiple_thetacs, but I
%broke it out so that we can run it multiple times with different l's,
%after loading the data

%% Confidence sets for linear combination of theta's
G_array = G_array_cell{1};
mean_g = mean(G_array(1,:,1));

%If l not specified, then do the l that gives you the mean weights
if( exist('l') == 0)
l = [ones(num_F_groups_parameters,1) / num_F_groups_parameters; mean_g];
%l = [1;0];
end

delta_true = [repmat(theta_c_true,num_F_groups_parameters,1);  theta_g_true];
l'*delta_true

%load( strcat( data_output_dir, dirname, 'Interacted_Moments/values_for_lp') )
%load( strcat( data_output_dir, dirname, 'Interacted_Moments/grid_cell') )

%For the LF methods, we save bounds of confidence sets, since we compute
%these exactly. These are intiated here
confidence_sets_using_c_alpha = NaN(numdatasets,2);
confidence_sets_using_c_lp_alpha = NaN(numdatasets,2);

%For the conditional/hybrid methods, we get rejection probabilities over a
%grid of values for l*theta. So we initialize those grids here
num_beta0_gridpoints = 201;
beta0_grid = linspace( xlim_graph(1), xlim_graph(2), num_beta0_gridpoints);
rejection_grid_conditional = NaN(numdatasets, num_beta0_gridpoints);
rejection_grid_hybrid = NaN(numdatasets, num_beta0_gridpoints);


nummoments = size( y_T_cell{ds,1} , 1);
Z_draws_interacted = randn(nummoments, 10000);

%%CHANGE THIS BACK
for ds = 1:numdatasets
%parfor ds = 1:numdatasets
    
   X_T = X_T_cell{ds,1};
   y_T = y_T_cell{ds,1};
   Sigma = Sigma_conditional_cell{ds,1}; 
   
   %Rescale moments by the variance (this is done earlier now)
%    D_sigma_minushalf = diag( sqrt( diag(Sigma) ).^(-1) );
%    y_T = D_sigma_minushalf * y_T;
%    X_T = D_sigma_minushalf * X_T;
%    
   
   %Do the LF test and store the resulting confidence set
   c_alpha = c_lf(Sigma, alpha, Z_draws_interacted); 
   confidence_sets_using_c_alpha(ds,:) = cs_linear_delta_lp_fn(y_T,X_T,l,c_alpha)';
   
   %Do the modified least favorable test and store the resulting confidence
   %set
   c_lp_alpha = c_lf_lp(X_T,Z_draws_interacted(:,1:numsims_lp),Sigma,alpha);
   confidence_sets_using_c_lp_alpha(ds,:) = cs_linear_delta_lp_fn(y_T,X_T,l,c_lp_alpha)';
   
   %%%Do the conditional and hybrid tests treating l as a non-linear
   %parameter
   
   %First, transform the moments slightly to avoid non-uniqueness
   %Transfrom the moments by adding in 1/100,000 times the average moment, so that solution to LP is unique
   c = 0.00001;
   [y_T,X_T,Sigma] = add_mean_to_moments_fn(c, y_T,X_T,Sigma);

   %We do a change of basis such that M*theta yields l*theta in the first
   %row (I construct an orthonormal basis for the rest of the basis
   %vectors, although this is not strictly necessary)
   M = [l' ; null( repmat(l',size(l)) )'];
   X_T = X_T * M^(-1);
   
   
    
    
    conditional_rejection_vec = NaN(num_beta0_gridpoints,1);
    hybrid_rejection_vec = NaN(num_beta0_gridpoints,1);
    
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
        
        
        
        count = count +1;
    end
    
    rejection_grid_conditional(ds,:) = conditional_rejection_vec;
    rejection_grid_hybrid(ds,:) = hybrid_rejection_vec;
    
end

%Create the confidence sets using the grid search
%grid_min_max
    
    %Average the rejection probabilities for the hybrid and conditional
    %methods
    rejection_grid_conditional = mean(rejection_grid_conditional,1);
    rejection_grid_hybrid = mean(rejection_grid_hybrid,1);
    

    ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/confidence_sets_lp');
    mkdir(ds_name);
    save( ds_name, 'confidence_sets_using_c_alpha', 'confidence_sets_using_c_lp_alpha',...
                   'rejection_grid_hybrid', 'rejection_grid_conditional','beta0_grid' );

    
    
    
%% Estimate the bounds of the identified set by taking the whole chain and setting critical value to 0


display('Starting to find identified set');

    long_ds_object = load( char(strcat( data_input_dir, dirname, 'ds_long.mat') )) ;
    
    
    F_array = long_ds_object.F_array;
    G_array = long_ds_object.G_array;
    Eta_jt_shocks_array = long_ds_object.Eta_jt_shocks_array;
    Eta_t_vec = long_ds_object.Eta_t_vec;
    %Pi_star_array = long_ds_object.Pi_star_array;
    J_t_array = long_ds_object.J_t_array;
    J_tminus1_array = long_ds_object.J_tminus1_array;
    Pi_array = long_ds_object.Pi_array;

    clear long_ds_object;
    
     [A_g_cell, A_c_cell, Y_cell] = generate_moment_fn_multiple_thetacs( F_group_cell_moments, F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array, use_basic_moments);

        first_iter = 1;
        parameter_number = 1;
        
      for(i = 1:num_F_groups_moments)
           
            %Update parameter number
            if( ~ismember( F_group_cell_moments{i},...
                           F_group_cell_parameters{parameter_number} ) )
                       parameter_number = parameter_number + 1;
            end
                
                
            %Get A_g and A_c as a function of lambda from the cell
            A_g_fn = A_g_cell{i,1};
            A_c_fn = A_c_cell{i,1};
            
            %Evaluate at lambda (default is true lambda)
            A_g = A_g_fn(lambda);
            A_c = A_c_fn(lambda);

            %Get Y from the cell
            Y = Y_cell{i,1};
            
            
            % To get the conditional matrix, we want to construct a matrix
            % A so that each row is all of the relevant weights on theta_g
            % and the theta_c's for a given market.
            % We construct the parts related to theta_c and theta_g
            % separately
            if(first_iter == 1)
                A_c_combined = A_c;
                A_g_combined = A_g;
                
            else
                A_c_combined = [A_c_combined, A_c];
                A_g_combined = [A_g_combined, A_g];
            end
            
            
            %We want to be able to write the moments as y_T - X_T * delta
            
            %Let Ybar_i be the sum (down columns) of Y for F-group i (that
            %is, summing the Y's for each moment over all of the markets)
            
            %Similarly, let X_i be the sum down columns for any variable X
            %from F-group i
            
            % We create two matrices:
           
            % Y = [ Ybar_1;...
            %      ;Ybar_NF]
            % X = [ Abar_c1, 0 ..., Abar_g1,
            %       0      , Abar_c2,...,Abar_g1 ...]
            
            A_c_bar = sum( A_c, 1)';
            A_g_bar = sum( A_g, 1)';
            Y_bar = sum(Y, 1)';
            
            length_A = size(A_c_bar,1);
            num_remaining_groups = num_F_groups_parameters - parameter_number;
            if(first_iter ==1)
                X_T = [ A_c_bar, zeros(length_A, num_remaining_groups), A_g_bar];
                
                y_T = -Y_bar;
                
                Y_wide = Y;
                
            else
                num_previous_groups = parameter_number -1;
                new_row_X = [zeros(length_A, num_previous_groups), A_c_bar, zeros(length_A, num_remaining_groups), A_g_bar];
                X_T = [X_T;new_row_X];
                
                y_T = [y_T;-Y_bar];
                
                Y_wide = [Y_wide, Y];
            end
            
            
            first_iter = 0;
        end
         
          
       
        %If combined_theta_g_moments, combined all of the theta_g moments
        %in to one set, rather than having one for each firm group
        
        if( combine_theta_g_moments == 1)        
            
            moment_nums = 1:size(X_T,1);
            moment_num_in_group = mod2( moment_nums, size(X_T,1) / num_F_groups_moments );
            theta_g_cols = moment_num_in_group == 5 | moment_num_in_group == 6;            
            moment_nums( theta_g_cols ) = moment_num_in_group( theta_g_cols);

            A_c_combined = grpstats2( A_c_combined' , moment_nums')';
            A_g_combined = grpstats2( A_g_combined' , moment_nums')';
            Y_wide = grpstats2( Y_wide', moment_nums')'; 
            X_T = grpstats2( X_T , moment_nums');
            y_T = grpstats2( y_T , moment_nums');
            
        end
        
        A = [A_c_combined , A_g_combined];
        %Remove any all zero columns
        A= A(:,any(A));
          
    
  %y_T and x_T are constructed so that population moments = y_T - X_T * delta
    % WE construct these so that they are less than 0 in expectation (the
    % moment fns are constructed so that y_T + X_T * delta is greater than 0 in expectation)
  T = size(A_g,1);
  X_T = X_T / sqrt( T ); 
  y_T = y_T / sqrt(T);
  


  %We set c_alpha = 0
  identified_set_bounds = cs_linear_delta_lp_fn(y_T,X_T,l,0)';

    ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/identified_set_bounds');
    mkdir(ds_name);
    save( ds_name, 'identified_set_bounds');

  display('Done with identified set calc');
