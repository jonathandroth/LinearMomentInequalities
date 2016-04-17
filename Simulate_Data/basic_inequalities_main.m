% This script does power calcs using the basic techniques for two
% parameters, hodling one fixed
tic;
alpha = 0.05;
beta = 0.005;


numgridpoints = 21;

lambda_grid = linspace(0,1,numgridpoints);
theta_c_grid = linspace(80,180,numgridpoints);
theta_g_grid = linspace(-75,25, numgridpoints);
% theta_c_grid = linspace(40,220,numgridpoints);
% theta_g_grid = linspace(-125,25, numgridpoints);


lambda_true = 0.386;
theta_c_true = 129.73;
theta_g_true = -21.38;



%Create a matrix of normal draws with mean 0 and covariace sigma
%Each column is a draw
nummoments = 6;

%Create a matrix of standard normals of size k x 10000 for simulating
%critical values
Z_draws = randn(nummoments, 10000);

numdatasets = 10;


for ds = 1:numdatasets
    
    full_moment_fn = generate_moment_fn( ds); 
    moment_fn = @(theta_c, theta_g) -full_moment_fn(theta_c, theta_g, lambda_true); 
    
    [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetac_thetag( moment_fn, theta_c_grid, theta_g_grid,Z_draws, alpha, beta);
 
    ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/grid', num2str(ds));
    save( ds_name, 'grid_lf', 'grid_rsw', 'grid_conditional', 'grid_hybrid');
end



%% Aggregate all the results to get the probabilities

rejection_prob_lf = zeros( numgridpoints);
rejection_prob_rsw = rejection_prob_lf;
rejection_prob_conditional = rejection_prob_lf;
rejection_prob_hybrid = rejection_prob_lf;

for ds = 1:numdatasets

    ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/grid', num2str(ds));
    load(ds_name);
    
    rejection_prob_lf = rejection_prob_lf + grid_lf / numdatasets;
    rejection_prob_rsw = rejection_prob_rsw + grid_rsw / numdatasets;
    rejection_prob_conditional = rejection_prob_conditional + grid_conditional / numdatasets;
    rejection_prob_hybrid = rejection_prob_hybrid + grid_hybrid / numdatasets;

end
[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);
 
contour( xgrid, ygrid, rejection_prob_conditional, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
contourf( xgrid, ygrid, rejection_prob_conditional);

toc;