
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
  