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
nummoments_basic = 6;
nummoments_interacted = 18;

%Create a matrix of standard normals of size k x 10000 for simulating
%critical values
Z_draws_interacted = randn(nummoments_interacted, 10000);
Z_draws_basic = Z_draws_interacted( 1:nummoments_basic,:);

numdatasets = 2;

dirnames = { 'Calibrated_SigmaZeta/', 'Calibrated_SigmaZeta_Over4/', 'SigmaZeta_Equal0/'};
    

for dirname = dirnames

    for ds = 1:numdatasets
    
    input_dir = strcat( '../../Output/Simulated_Data/', dirname)
    ds
    
    [moment_fn_allparams,moment_fn_interacted_allparams] = generate_moment_fn( input_dir, ds); 
    
    %Compute and save the grids for the basic moments
    moment_fn = @(theta_c, theta_g) -moment_fn_allparams(theta_c, theta_g, lambda_true);
    [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetac_thetag( moment_fn, theta_c_grid, theta_g_grid,Z_draws_basic, alpha, beta);
 
    ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, 'Basic_Moments/grid', num2str(ds));
    ds_name = ds_name{:};
    save( ds_name, 'grid_lf', 'grid_rsw', 'grid_conditional', 'grid_hybrid');
    
    %Compute and save the grid for the interacted moments
    moment_fn_interacted = @(theta_c, theta_g) -moment_fn_interacted_allparams(theta_c, theta_g, lambda_true);
    [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetac_thetag( moment_fn_interacted, theta_c_grid, theta_g_grid,Z_draws_interacted, alpha, beta);
    
    ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, 'Interacted_Moments/grid', num2str(ds));
    ds_name = ds_name{:};
    save( ds_name, 'grid_lf', 'grid_rsw', 'grid_conditional', 'grid_hybrid');
    end

end

%% Aggregate all the results to get the probabilities
for moment_type = {'Basic_Moments/', 'Interacted_Moments/'};
    
for dirname = dirnames
    
rejection_prob_lf = zeros( numgridpoints);
rejection_prob_rsw = rejection_prob_lf;
rejection_prob_conditional = rejection_prob_lf;
rejection_prob_hybrid = rejection_prob_lf;

for ds = 1:numdatasets

    ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, moment_type, 'grid', num2str(ds));
    ds_name = ds_name{:};
    load(ds_name);
    
    rejection_prob_lf = rejection_prob_lf + grid_lf / numdatasets;
    rejection_prob_rsw = rejection_prob_rsw + grid_rsw / numdatasets;
    rejection_prob_conditional = rejection_prob_conditional + grid_conditional / numdatasets;
    rejection_prob_hybrid = rejection_prob_hybrid + grid_hybrid / numdatasets;

end
[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);

figure_dir = strcat( '../../Output/Figures/Rejection_Grids/Lambda_Constant/', dirname, moment_type);
figure_dir = figure_dir{:};


contour( xgrid, ygrid, rejection_prob_lf, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
title('Rejection Probabilities - Least Favorable');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_dir, 'rejection_probs_lf'), 'epsc');

contour( xgrid, ygrid, rejection_prob_rsw, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
title('Rejection Probabilities - RSW');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_dir, 'rejection_probs_rsw'), 'epsc');

contour( xgrid, ygrid, rejection_prob_hybrid, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
title('Rejection Probabilities - Hybrid');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_dir, 'rejection_probs_hybrid'), 'epsc');

contour( xgrid, ygrid, rejection_prob_conditional, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
title('Rejection Probabilities - Conditional');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_dir, 'rejection_probs_conditional'), 'epsc');


% contourf( xgrid, ygrid, rejection_prob_conditional);
% contour( xgrid, ygrid, rejection_prob_conditional, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');

end
end
toc;