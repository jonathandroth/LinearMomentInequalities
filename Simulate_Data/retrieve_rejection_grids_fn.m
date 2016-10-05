function [rejection_prob_lf, rejection_prob_rsw, rejection_prob_conditional, rejection_prob_hybrid] = ...
    retrieve_rejection_grids_fn( rejection_grids_cell,theta_c_grid, theta_g_grid )

numdatasets = size( rejection_grids_cell, 1);

rejection_prob_lf = zeros( length(theta_c_grid), length(theta_g_grid));
rejection_prob_rsw = rejection_prob_lf;
rejection_prob_conditional = rejection_prob_lf;
rejection_prob_hybrid = rejection_prob_lf;

for ds = 1:numdatasets


    grid_lf = rejection_grids_cell{ds,1};
    grid_rsw = rejection_grids_cell{ds,2};
    grid_conditional = rejection_grids_cell{ds,3};
    grid_hybrid = rejection_grids_cell{ds,4};
    
    
    rejection_prob_lf = rejection_prob_lf + grid_lf / numdatasets;
    rejection_prob_rsw = rejection_prob_rsw + grid_rsw / numdatasets;
    rejection_prob_conditional = rejection_prob_conditional + grid_conditional / numdatasets;
    rejection_prob_hybrid = rejection_prob_hybrid + grid_hybrid / numdatasets;
    
    
end