%% Specify parameters
lp_set_parameters
%% Import data and calculate moments and covariance matrices
%lp_create_moments_and_covariances


lambda_vec = (0.01:.2:5.01)';
%lambda_vec = (0.01:5:5.01)';

conditional_test_noadjustment = NaN(numdatasets, size(lambda_vec,1));
conditional_test_adjustment = NaN(size(conditional_test_noadjustment));
hybrid_test_adjustment = NaN(size(conditional_test_noadjustment));
lf_test_original = NaN(size(conditional_test_noadjustment));
lf_test_modified = NaN(size(conditional_test_noadjustment));

identified_set = NaN( size(lambda_vec) );

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


conditional_test_noadjustment(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);

%Transfrom the moments by adding in 1/100,000 times the average moment, so that solution to LP is unique
c = 0.00001;
k = size(X_T,1);
A = ( eye(k,k) + c*ones(k,k) );

X_T = A * X_T;
y_T = A * y_T;
Sigma = A * Sigma * A';

%Renormalize the moments
D_minushalf = diag( diag(Sigma).^(-1/2) );
X_T = D_minushalf * X_T;
y_T = D_minushalf * y_T;
Sigma = D_minushalf * Sigma * D_minushalf;

conditional_test_adjustment(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);
hybrid_test_adjustment(ds,i) = lp_hybrid_test_fn( y_T, X_T, Sigma, alpha, alpha/10);

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
    save( strcat(ds_dir, 'lambda_results'), 'conditional_test_adjustment', 'conditional_test_noadjustment', 'lf_test_original', 'lf_test_modified', 'hybrid_test_adjustment');


 %% Estimate the identified set for lambda using a large chain and seeing where the moments are <= 0
    display('Starting to find identified set for lambda');

    
   
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

lambda_index = 1;     
 for lambda = lambda_vec'
        lambda
        
        first_iter = 1;
        for(i = 1:num_F_groups)
           
            %Get A_g and A_c as a function of lambda from the cell
            A_g_fn = A_g_cell{i,1};
            A_c_fn = A_c_cell{i,1};
            
            %Evaluate at the current lambda
            A_g = A_g_fn(lambda);
            A_c = A_c_fn(lambda);

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
            num_remaining_groups = num_F_groups - i;
            if(first_iter ==1)
                X_T = [ A_c_bar, zeros(length_A, num_remaining_groups), A_g_bar];
                
                y_T = -Y_bar;
                
                Y_wide = Y;
                
            else
                num_previous_groups = i -1;
                new_row_X = [zeros(length_A, num_previous_groups), A_c_bar, zeros(length_A, num_remaining_groups), A_g_bar];
                X_T = [X_T;new_row_X];
                
                y_T = [y_T;-Y_bar];
                
                Y_wide = [Y_wide, Y];
            end
            
            
            first_iter = 0;
        end
         
          
        if( combine_theta_g_moments == 1)        
            moment_nums = 1:size(A,2);
            moment_nums_mod6 = mod2( moment_nums,6 );
            theta_g_cols = moment_nums_mod6 == 5 | moment_nums_mod6 == 6;       
            moment_nums_theta_g_cols = mod2( moment_nums( theta_g_cols), size(A,2) / num_F_groups);
            moment_nums( theta_g_cols) = moment_nums_theta_g_cols;

            A = grpstats2( A' , moment_nums')';

            moment_nums = 1:size(X_T,1);
            moment_nums_mod6 = mod2( moment_nums,6 );
            theta_g_cols = moment_nums_mod6 == 5 | moment_nums_mod6 == 6;       
            moment_nums_theta_g_cols = mod2( moment_nums( theta_g_cols), size(X_T,1) / num_F_groups);
            moment_nums( theta_g_cols) = moment_nums_theta_g_cols;

            X_T = grpstats2( X_T , moment_nums');
            y_T = grpstats2( y_T , moment_nums');
            Y_wide = grpstats2( Y_wide', moment_nums')'; 
 
        end
   
        
        
    
  %y_T and x_T are constructed so that population moments = y_T - X_T * delta
    % WE construct these so that they are less than 0 in expectation (the
    % moment fns are constructed so that y_T + X_T * delta is greater than 0 in expectation)
  T = size(A_g,1);
  X_T = X_T / sqrt( T ); 
  y_T = y_T / sqrt(T);
  
    %Transform the moments to add 10^(-5) times the mean
    c = 0.00001;    
    k = size(X_T,1);
    A = ( eye(k,k) + c*ones(k,k) );

    X_T = A * X_T;
    y_T = A * y_T;

  
  %We say lambda is in the identified set if there is any delta such that
  %the moments hold in our long chain, i.e. if eta <= 0
  eta  = test_delta_lp_fn( y_T, X_T, optimoptions('linprog','Algorithm','dual-simplex', 'Display', 'off'));
  identified_set(lambda_index) = (eta <= 0 );
  
  lambda_index = lambda_index+1;
 end
 
 
 save( strcat(ds_dir, 'lambda_identified_set'), 'identified_set', 'lambda_vec');

 %% Plot

 ds_dir = strcat( data_output_dir, 'Interacted_Moments/');
 load(  strcat(ds_dir, 'lambda_identified_set') );
 load(  strcat(ds_dir, 'lambda_results') );
 
  
 plot(lambda_vec, [mean(conditional_test_adjustment);...
                  mean(conditional_test_noadjustment);...
                  mean(hybrid_test_adjustment);...
                  mean(lf_test_original);...
                  mean(lf_test_modified)] ) 

identified_set_max = max( lambda_vec( identified_set == 1) );

line( [identified_set_max; identified_set_max], [0;1], 'LineStyle', '--', 'Color',  'r');

legend( 'Conditional (Transformed)', 'Conditional (Non-Transformed)', 'Hybrid', 'LF','LF (modified)', 'Identified Set Bound', 'Location','eastoutside' );
ylabel('Rejection Probability');
xlabel('Lambda');
title('Rejection Probabilities for Lambda');
saveas( gcf, strcat(figures_output_dir,filename_graph ), 'epsc');

