% This script does power calcs using the basic techniques for two
% parameters, hodling one fixed
tic;
alpha = 0.05;
beta = 0.005;


numgridpoints = 21;

lambda_grid = linspace(0,1,numgridpoints);
theta_c_grid = linspace(80,180,numgridpoints);
theta_g_grid = linspace(-75,25, numgridpoints);

lambda_true = 0.386;
theta_c_true = 129.73;
theta_g_true = -21.38;


rejection_prob_lf = zeros( numgridpoints);
rejection_prob_rsw = rejection_prob_lf;
rejection_prob_conditional = rejection_prob_lf;
rejection_prob_hybrid = rejection_prob_lf;

%Create a matrix of normal draws with mean 0 and covariace sigma
%Each column is a draw
nummoments = 6;

%Create a matrix of standard normals of size k x 10000 for simulating
%critical values
Z_draws = randn(nummoments, 10000);

numdatasets = 50;


for ds = 1:numdatasets
    
    full_moment_fn = generate_moment_fn( ds); 
    moment_fn = @(theta_c, theta_g) -full_moment_fn(theta_c, theta_g, lambda_true); 
    
for theta_c_index = 1:numgridpoints
    for theta_g_index = 1:numgridpoints
        
        theta_c = theta_c_grid(theta_c_index);
        theta_g = theta_g_grid(theta_g_index);
        
        moments_mat = moment_fn(theta_c, theta_g);
        g_T = 1/ sqrt( size(moments_mat,1) ) * sum(moments_mat,1)';
        Sigma = cov(moments_mat);
        
        %Do the tests
        [test_lf, test_rsw, test_conditional, test_hybrid] = basic_tests(g_T, Sigma, Z_draws, alpha, beta);
        
        %Update the rejection probability matrices
        rejection_prob_lf(theta_c_index, theta_g_index) =  rejection_prob_lf(theta_c_index, theta_g_index) + test_lf;
        rejection_prob_rsw(theta_c_index, theta_g_index) =  rejection_prob_rsw(theta_c_index, theta_g_index) + test_rsw;
        rejection_prob_conditional(theta_c_index, theta_g_index) =  rejection_prob_conditional(theta_c_index, theta_g_index) + test_conditional;
        rejection_prob_hybrid(theta_c_index, theta_g_index) =  rejection_prob_hybrid(theta_c_index, theta_g_index) + test_hybrid;
    end
    
end

end

rejection_prob_lf = rejection_prob_lf / numdatasets;
rejection_prob_rsw = rejection_prob_rsw / numdatasets;
rejection_prob_conditional = rejection_prob_conditional / numdatasets;
rejection_prob_hybrid = rejection_prob_hybrid / numdatasets;

[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);
 
contour( xgrid, ygrid, rejection_prob_conditional, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
contourf( xgrid, ygrid, rejection_prob_conditional);

toc;