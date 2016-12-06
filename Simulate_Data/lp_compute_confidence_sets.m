%This script computes the confidence sets and identified set for l * theta
%in the lp_confidence_sets script

%If l is not specfied in the environment, it is assumed to be [1,0,...,0,
%mean_g]

%This script used to be part of lp_confidence_sets_multiple_thetacs, but I
%broke it out so that we can run it multiple times with different l's,
%after loading the data

%% Confidence sets for linear combination of theta's
G_array = G_array_cell{1};
mean_g = mean(G_array(1,:,1));

%If l not specified, then do the l that gives you the mean weights
if( exist('l') == 0)
l = [1; zeros(num_F_groups-1,1); mean_g];
%l = [1;0];
end

delta_true = [repmat(theta_c_true,num_F_groups,1);  theta_g_true];
l'*delta_true

load( strcat( data_output_dir, dirname, 'Interacted_Moments/values_for_lp') )
%load( strcat( data_output_dir, dirname, 'Interacted_Moments/grid_cell') )

confidence_sets_using_c_alpha = NaN(numdatasets,2);
confidence_sets_using_c_lp_alpha = NaN(numdatasets,2);

nummoments = size( y_T_cell{ds,1} , 1);
Z_draws_interacted = randn(nummoments, 10000);

parfor ds = 1:numdatasets
    
   X_T = X_T_cell{ds,1};
   y_T = y_T_cell{ds,1};
   Sigma = Sigma_conditional_cell{ds,1}; 
   
   D_sigma_minushalf = diag( sqrt( diag(Sigma) ).^(-1) );
   y_T = D_sigma_minushalf * y_T;
   X_T = D_sigma_minushalf * X_T;
   
   
   c_alpha = c_lf(Sigma, alpha, Z_draws_interacted);
   
   %test_delta_lp_fn(y_T, X_T);
   
 
   confidence_sets_using_c_alpha(ds,:) = cs_linear_delta_lp_fn(y_T,X_T,l,c_alpha)';
   
   c_lp_alpha = c_lf_lp(X_T,Z_draws_interacted(:,1:numsims_lp),Sigma,alpha);
   confidence_sets_using_c_lp_alpha(ds,:) = cs_linear_delta_lp_fn(y_T,X_T,l,c_lp_alpha)';
end

%Create the confidence sets using the grid search
%grid_min_max


    ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/confidence_sets_lp');
    mkdir(ds_name);
    save( ds_name, 'confidence_sets_using_c_alpha', 'confidence_sets_using_c_lp_alpha');

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
    
     [A_g_cell, A_c_cell, Y_cell] = generate_moment_fn_multiple_thetacs( F_group_cell, F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array);

        first_iter = 1;
        for(i = 1:num_F_groups)
           
            %Get A_g and A_c as a function of lambda from the cell
            A_g_fn = A_g_cell{i,1};
            A_c_fn = A_c_cell{i,1};
            
            %Evaluate at the true lambda
            A_g = A_g_fn(lambda_true);
            A_c = A_c_fn(lambda_true);

            %Get Y from the cell
            Y = Y_cell{i,1};
            
            
            % To get the conditional matrix, we want to construct a matrix
            % A so that each row is all of the relevant weights on theta_g
            % and the theta_c's for a given market
            if(first_iter == 1)
                A = [A_c , A_g];
                
            else
                A = [A, A_c, A_g];
            end
            
            
            %We want to be able to write the moments as Y_T - X_T * delta
            
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
            num_remaining_groups = num_F_groups - i;
            if(first_iter ==1)
                X_T = [ A_c_bar, zeros(length_A, num_remaining_groups), A_g_bar];
                
                Y_T = -Y_bar;
                
                Y_wide = Y;
                
            else
                num_previous_groups = i -1;
                new_row_X = [zeros(length_A, num_previous_groups), A_c_bar, zeros(length_A, num_remaining_groups), A_g_bar];
                X_T = [X_T;new_row_X];
                
                Y_T = [Y_T;-Y_bar];
                
                Y_wide = [Y_wide, Y];
            end
            
            
            first_iter = 0;
        end
         
          
         
    
  %y_T and x_T are constructed so that population moments = y_T - X_T * delta
    % WE construct these so that they are less than 0 in expectation (the
    % moment fns are constructed so that y_T + X_T * delta is greater than 0 in expectation)
  T = size(A_g,1);
  X_T = X_T / sqrt( T ); 
  y_T = Y_T / sqrt(T);
  
  %We set c_alpha = 0
  identified_set_bounds = cs_linear_delta_lp_fn(y_T,X_T,l,0)';

    ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/identified_set_bounds');
    mkdir(ds_name);
    save( ds_name, 'identified_set_bounds');

  display('Done with identified set calc');
