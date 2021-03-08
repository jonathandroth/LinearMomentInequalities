    %%% This script saves the Y component of the moments from the examples
    %%% with linear parameters. 
    
    %%
    use_basic_moments = 1;

%specify whether to combine the moments for theta_g, or to have separate
%ones for each firm
combine_theta_g_moments = 1;

%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/1Thetac/Basic_Moments/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/1Thetac/Basic_Moments/';


%Specify the groups of firms that have different coefficients
F_group_cell_moments = {[1;2;3;4;5;6;7;8;9]};


%num_F_groups_moments = size(F_group_cell_moments,1);
%l = [1; zeros(num_F_groups_moments-1,1); mean_g];

%Run the main script (this does confidence sets for the mean weight)
if( exist('xlim_graph_meanweight'))
    xlim_graph = xlim_graph_meanweight;
else
    xlim_graph = [83;84];
end

if( exist('xsplit_graph_meanweight'))
   xsplit_graph = xsplit_graph_meanweight;
end

filename_graph = 'Mean_Weight_Rejection_Probabilities';
lp_set_parameters;
    
  
    %%
    long_ds_object = load( char(strcat( data_input_dir, dirname, 'ds_long.mat') )) ;
    length_long_chain = size( long_ds_object.J_t_array,3);

    rand_index_oracle = randsample( length_long_chain , 10*nummarkets);
    
    
    F_array = long_ds_object.F_array(:,:,rand_index_oracle);
    G_array = long_ds_object.G_array(:,:,rand_index_oracle);
    Eta_jt_shocks_array = long_ds_object.Eta_jt_shocks_array(:,:,rand_index_oracle);
    Eta_t_vec = long_ds_object.Eta_t_vec(rand_index_oracle);
    %Pi_star_array = long_ds_object.Pi_star_array(:,:,rand_index);
    J_t_array = long_ds_object.J_t_array(:,:,rand_index_oracle);
    J_tminus1_array = long_ds_object.J_tminus1_array(:,:,rand_index_oracle);
    Pi_array = long_ds_object.Pi_array(:,:,rand_index_oracle);
    
   clear long_ds_object    
    
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
        
        Y_mat = - Y_wide;
        Y_star_wide = -Y_wide - theta_c_true * A_c_combined - theta_g_true * A_g_combined;
        full_dir_name = strcat( data_output_dir, dirname);
        mkdir(full_dir_name);
        ds_name = strcat(full_dir_name, 'y_i_values');
        save( ds_name, 'Y_star_wide', 'Y_mat', 'A_c_combined', 'A_g_combined');
    