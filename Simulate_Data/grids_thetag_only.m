%This function takes a moment_fn (a function of theta_g that generates the
%moments), a grid for theta_g, and Z_draws, alpha, beta.

%It returns a series of binary grids as to whether the given value of
%theta_g was rejected or accepted by the various tests

function [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetag_only( moment_fn, theta_g_grid,Z_draws, alpha, beta)

 grid_lf = NaN( size(theta_g_grid) );
 grid_rsw = grid_lf;
 grid_conditional = grid_lf;
 grid_hybrid = grid_lf;
 
   
    for theta_g_index = 1:length(theta_g_grid)
        
      
        theta_g = theta_g_grid(theta_g_index);
        
        moments_mat = moment_fn(theta_g);
        g_T = 1/ sqrt( size(moments_mat,1) ) * sum(moments_mat,1)';
        Sigma = cov(moments_mat);
        
        %Do the tests
        [test_lf, test_rsw, test_conditional, test_hybrid] = basic_tests(g_T, Sigma, Z_draws, alpha, beta);
        
        %Update the rejection probability matrices
        grid_lf(theta_g_index) =  test_lf;
        grid_rsw(theta_g_index) =  test_rsw;
        grid_conditional(theta_g_index) =  test_conditional;
        grid_hybrid(theta_g_index) =  test_hybrid;
    end   
end
