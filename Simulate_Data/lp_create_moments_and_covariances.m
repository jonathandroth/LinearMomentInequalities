
%%This script imports the datasets and computes and stores the moment functions and
%%conditional variance matrices for each of them

%This script is called during the beginning of the lp_confidence_sets_script, after
%specifying the basic parameters. I've broken this into a separate script
%so that it can be run only once, and then the confidence sets can be
%calculated for multiple linear combinations afterwards


%% Calculate the moments and conditional variances

timing_vec_lf = zeros(numdatasets,1);
timing_vec_lfp = zeros(numdatasets,1);
timing_vec_conditional = zeros(numdatasets,1);
timing_vec_hybrid = zeros(numdatasets,1);
timing_vec_as = zeros(numdatasets,1);
timing_vec_kms = zeros(numdatasets,1);
timing_vec_rcc = zeros(numdatasets,1);
timing_vec_cc = zeros(numdatasets,1);


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

    rng(ds + dsoffset);
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
   
    %ds
    
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
        
        startTimeCondVar = tic;
        %Calculate the variance of Y conditional on A using the Abadie et
        %al matched pairs method
        Sigma_conditional = conditional_variance_fn(Y_wide, A, diagonal);
        runTimeCondVar = toc(startTimeCondVar);
        
        
        
        %Add time for conditional variance to all the timers
        timing_vec_lf(ds) = timing_vec_lf(ds) + runTimeCondVar;
        timing_vec_lfp(ds) = timing_vec_lfp(ds) + runTimeCondVar;
        timing_vec_conditional(ds) = timing_vec_conditional(ds) + runTimeCondVar;
        timing_vec_hybrid(ds) = timing_vec_hybrid(ds) + runTimeCondVar;
        timing_vec_as(ds) = timing_vec_as(ds) + runTimeCondVar;
        timing_vec_kms(ds) = timing_vec_kms(ds) + runTimeCondVar;
        timing_vec_rcc(ds) = timing_vec_rcc(ds) + runTimeCondVar;
        timing_vec_cc(ds) = timing_vec_cc(ds) + runTimeCondVar;

        
  %y_T and x_T are constructed so that population moments = y_T - X_T * delta
    % WE construct these so that they are less than 0 in expectation (the
    % moment fns are constructed so that y_T + X_T * delta is greater than 0 in expectation)
  T = size(A_g,1);
  X_T = X_T / sqrt( T ); 
  y_T = y_T / sqrt(T);
  %max( abs( (y_T + X_T * [1 ;2]) - mean(  moment_fn_interacted_allparams(1,2, lambda) )' ))
  
  %Normalize the moments and covariance matrix so that Sigma has diagnol of
  %1
  D_sigma_minushalf = diag( sqrt( diag(Sigma_conditional) ).^(-1) );
  y_T = D_sigma_minushalf * y_T;
  X_T = D_sigma_minushalf * X_T;
  Sigma_conditional = D_sigma_minushalf * Sigma_conditional * D_sigma_minushalf';
  
  
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
    mkdir(ds_name);
    save( ds_name, 'y_T_cell', 'X_T_cell', 'Sigma_conditional_cell', 'c_lf_cell');
    


toc;
