function [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetac_thetag( moment_fn, theta_c_grid, theta_g_grid,Z_draws, alpha, beta)

 grid_lf = NaN( length(theta_c_grid), length(theta_g_grid) );
 grid_rsw = grid_lf;
 grid_conditional = grid_lf;
 grid_hybrid = grid_lf;
   
for theta_c_index = 1:length(theta_c_grid)
    for theta_g_index = 1:length(theta_g_grid)
        
        theta_c = theta_c_grid(theta_c_index);
        theta_g = theta_g_grid(theta_g_index);
        
        moments_mat = moment_fn(theta_c, theta_g);
        g_T = 1/ sqrt( size(moments_mat,1) ) * sum(moments_mat,1)';
        Sigma = cov(moments_mat);
        
        %Do the tests
        [test_lf, test_rsw, test_conditional, test_hybrid] = basic_tests(g_T, Sigma, Z_draws, alpha, beta);
        
        %Update the rejection probability matrices
        grid_lf(theta_c_index, theta_g_index) =  test_lf;
        grid_rsw(theta_c_index, theta_g_index) =  test_rsw;
        grid_conditional(theta_c_index, theta_g_index) =  test_conditional;
        grid_hybrid(theta_c_index, theta_g_index) =  test_hybrid;
    end   
end
