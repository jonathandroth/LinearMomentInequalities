function [A_g_cell, A_c_cell, Y_cell] = generate_moment_fn_multiple_thetacs( F_group_cell_moments, F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array, varargin)

%If an optional last argument is specified, it is an indicator for whether
%to use the basic moments only. Otherwise, assume this is zero and use the
%full interacted set of moments
if(isempty( varargin) == 0)
    use_basic_moments = varargin{1};
else
    use_basic_moments = 0;
end

total_num_firms = size(F_array,1);     %Use the fact that rows index firms, so that the first 
num_f_groups = size(F_group_cell_moments,1);



A_g_cell = cell(num_f_groups,1);
A_c_cell = cell(num_f_groups,1);
Y_cell = cell(num_f_groups,1);

for(f_group_index = 1:num_f_groups)

    f_group = F_group_cell_moments{f_group_index, 1};
    

    subset_indicator = ismember(1:total_num_firms, f_group);
    
    

[~,~,Y, A_g, A_c,Y_basic, A_g_basic, A_c_basic] = ...
        generate_moment_fn( F_array(subset_indicator,:,:), G_array(subset_indicator,:,:),...
        Eta_jt_shocks_array(subset_indicator,:,:), Eta_t_vec, Pi_array(subset_indicator, : ,:), ...
        J_t_array(subset_indicator, : ,:), J_tminus1_array(subset_indicator,:,:)); 
        
    
    if( use_basic_moments == 1)
        A_g_cell{f_group_index, 1} = A_g_basic;
        A_c_cell{f_group_index, 1} = A_c_basic;
        Y_cell{f_group_index, 1} = Y_basic;
        
    else
        A_g_cell{f_group_index, 1} = A_g;
        A_c_cell{f_group_index, 1} = A_c;
        Y_cell{f_group_index, 1} = Y;
    
    end
    
end
