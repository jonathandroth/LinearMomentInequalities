function [A_g_cell, A_c_cell, Y_cell] = generate_moment_fn_multiple_thetacs( F_group_cell, F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array)

total_num_firms = size(F_array,1);     %Use the fact that rows index firms, so that the first 
num_f_groups = size(F_group_cell,1);



A_g_cell = cell(num_f_groups,1);
A_c_cell = cell(num_f_groups,1);
Y_cell = cell(num_f_groups,1);

for(f_group_index = 1:num_f_groups)

    f_group = F_group_cell{f_group_index, 1};
    

    subset_indicator = ismember(1:total_num_firms, f_group);
    
    
    
[moment_fn_allparams,moment_fn_interacted_allparams,Y, A_g, A_c,Y_basic, A_g_basic, A_c_basic] = ...
        generate_moment_fn( F_array(subset_indicator,:,:), G_array(subset_indicator,:,:),...
        Eta_jt_shocks_array(subset_indicator,:,:), Eta_t_vec, Pi_array(subset_indicator, : ,:), ...
        J_t_array(subset_indicator, : ,:), J_tminus1_array(subset_indicator,:,:)); 

    A_g_cell{f_group_index, 1} = A_g;
    A_c_cell{f_group_index, 1} = A_c;
    Y_cell{f_group_index, 1} = Y;
    
    
    
end