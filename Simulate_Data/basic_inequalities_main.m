% This script does power calcs using the basic techniques for two
% parameters, hodling one fixed



%pc = parcluster('local');
%pc.NumWorkers = 2;

tic;
alpha = 0.05;
beta = 0.005;


numgridpoints = 10;

lambda_grid = linspace(0,1,numgridpoints);
%theta_c_grid = linspace(80,180,numgridpoints);
%theta_g_grid = linspace(-100,50, numgridpoints);
% theta_c_grid = linspace(40,220,numgridpoints);
% theta_g_grid = linspace(-125,25, numgridpoints);
theta_g_grid = -150:5:100;
theta_c_grid = -250:10:510;

lambda_true = 0.386;
theta_c_true = 129.73;
theta_g_true = -21.38;



%Create a matrix of normal draws with mean 0 and covariace sigma
%Each column is a draw
nummoments_basic = 6;
nummoments_interacted = 18;

%Create a matrix of standard normals of size k x 10000 for simulating
%critical values
rng(0);
Z_draws_interacted = randn(nummoments_interacted, 10000);
Z_draws_basic = Z_draws_interacted( 1:nummoments_basic,:);

numdatasets = 500;

nummarkets = 27; %This is the number of markets to sample from the long chain


%dirnames = { 'Calibrated_SigmaZeta/', 'Calibrated_SigmaZeta_Over4/', 'SigmaZeta_Equal0/'};
%dirnames = { 'Calibrated_SigmaZeta/'};
    
%parpool('local', 2);
    



for dirname = dirnames
    
    
    long_ds_object = load( char(strcat( data_input_dir, dirname, 'ds_long.mat') )) ;
    length_long_chain = size( long_ds_object.J_t_array,3);
    
    
    
    
    %Create cells to store the datasets and the rejection grids
    rejection_grids_cell = cell(numdatasets, 4) ;
    interacted_rejection_grids_cell = cell(numdatasets, 4) ;
    
    F_array_cell = cell(numdatasets,1);
    G_array_cell = cell(numdatasets,1);
    Eta_shocks_array_cell = cell(numdatasets,1);
    J_t_array_cell = cell(numdatasets,1);
    J_tminus1_array_cell = cell(numdatasets,1);
    Pi_array_cell = cell(numdatasets,1);

    %The following loop takes subsets of the long chain and stores in a
    %cell. Each individual dataset will then be passed to a worker in the
    %parfor loop below
    for ds = 1:numdatasets

    rng(ds);
    rand_index = randsample( length_long_chain , nummarkets);
    
    F_array_cell{ds} = long_ds_object.F_array(:,:,rand_index);
    G_array_cell{ds} = long_ds_object.G_array(:,:,rand_index);
    Eta_shocks_array_cell{ds} = long_ds_object.Eta_shocks_array(:,:,rand_index);
    %Pi_star_array = long_ds_object.Pi_star_array(:,:,rand_index);
    J_t_array_cell{ds} = long_ds_object.J_t_array(:,:,rand_index);
    J_tminus1_array_cell{ds} = long_ds_object.J_tminus1_array(:,:,rand_index);
    Pi_array_cell{ds} = long_ds_object.Pi_array(:,:,rand_index);
    end
    
    
    
    %If we want the "oracle covariance", then we take one long subset and
    %calculate the covariance (or conditional covariance) on that. 
    if( oracle_cov == 1)
        
        rng(0);
        rand_index = randsample( length_long_chain , 1000); %take a subset of length 1000
    
        F_array = long_ds_object.F_array(:,:,rand_index);
        G_array = long_ds_object.G_array(:,:,rand_index);
        Eta_shocks_array = long_ds_object.Eta_shocks_array(:,:,rand_index);
        %Pi_star_array = long_ds_object.Pi_star_array(:,:,rand_index);
        J_t_array = long_ds_object.J_t_array(:,:,rand_index);
        J_tminus1_array = long_ds_object.J_tminus1_array(:,:,rand_index);
        Pi_array = long_ds_object.Pi_array(:,:,rand_index);
    
         [moment_fn_allparams,moment_fn_interacted_allparams,Y, A_g, A_c,Y_basic, A_g_basic, A_c_basic] = ...
        generate_moment_fn( F_array, G_array, Eta_shocks_array, Pi_array, J_t_array, J_tminus1_array); 
        moment_fn = @(theta_c, theta_g) -moment_fn_allparams(theta_c, theta_g, lambda_true);

   
     if(conditional_cov == 1)
        %A_g and A_c come out as functions of labmda
        %We replace these with their value at the true lambda
        A_g = A_g(lambda_true);
        A_c = A_c(lambda_true); 
        A_g_basic = A_g_basic(lambda_true);
        A_c_basic = A_c_basic(lambda_true); 
        
        
        %Merge theta_g and theta_c coefficients
        A = [A_c, A_g];
        A_basic = [A_c_basic, A_g_basic];
        
        %Remove any all zero columns
        A= A(:,any(A));
        A_basic = A_basic(:,any(A_basic));
        
        %Calculate the variance of Y conditional on A using the Abadie et
        %al matched pairs method
        Sigma_conditional_oracle = conditional_variance_fn(Y, A, diagonal);
        Sigma_conditional_basic_oracle  =conditional_variance_fn(Y_basic, A_basic, diagonal);
     
 
     end
      
        %You can't do the unconditional here, because it depends on the
        %theta's. It must be done each time in the grid, if you want to do
        %it
     end
        
        
    
    
    clear long_ds_object    
    parfor ds = 1:numdatasets
        
    F_array = F_array_cell{ds};
    G_array = G_array_cell{ds};
    Eta_shocks_array = Eta_shocks_array_cell{ds};
    %Pi_star_array = Pi_star_array_cell{ds};
    J_t_array = J_t_array_cell{ds};
    J_tminus1_array = J_tminus1_array_cell{ds};
    Pi_array = Pi_array_cell{ds};
   
    ds
    
    [moment_fn_allparams,moment_fn_interacted_allparams,Y, A_g, A_c,Y_basic, A_g_basic, A_c_basic] = ...
        generate_moment_fn( F_array, G_array, Eta_shocks_array, Pi_array, J_t_array, J_tminus1_array); 
    moment_fn = @(theta_c, theta_g) -moment_fn_allparams(theta_c, theta_g, lambda_true);

   
    if(conditional_cov == 1) 
        if(oracle_cov == 0)
        %A_g and A_c come out as functions of labmda
        %We replace these with their value at the true lambda
        A_g = A_g(lambda_true);
        A_c = A_c(lambda_true); 
        A_g_basic = A_g_basic(lambda_true);
        A_c_basic = A_c_basic(lambda_true); 
        
        
        %Merge theta_g and theta_c coefficients
        A = [A_c, A_g];
        A_basic = [A_c_basic, A_g_basic];
        
        %Remove any all zero columns
        A= A(:,any(A));
        A_basic = A_basic(:,any(A_basic));
        
        %Calculate the variance of Y conditional on A using the Abadie et
        %al matched pairs method
        Sigma_conditional = conditional_variance_fn(Y, A, diagonal);
        Sigma_conditional_basic  =conditional_variance_fn(Y_basic, A_basic, diagonal);
         
        elseif(oracle_cov ==1)
            Sigma_conditional = Sigma_conditional_oracle;
            Sigma_conditional_basic = Sigma_conditional_basic_oracle;
        end
    end
    
    %Compute and save the grids for the basic moments
    if(conditional_cov == 1)
        [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetac_thetag( moment_fn, theta_c_grid, theta_g_grid,Z_draws_basic, alpha, beta, Sigma_conditional_basic);
    else
        [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetac_thetag( moment_fn, theta_c_grid, theta_g_grid,Z_draws_basic, alpha, beta);
    end
       
    rejection_grids_cell(ds,:) = {grid_lf, grid_rsw, grid_conditional,grid_hybrid};
    
%     ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, 'Basic_Moments/grid', num2str(ds));
%     ds_name = ds_name{:};
%     save( ds_name, 'grid_lf', 'grid_rsw', 'grid_conditional', 'grid_hybrid');
%     
    %Compute and save the grid for the interacted moments
    if(conditional_cov == 1)
        moment_fn_interacted = @(theta_c, theta_g) -moment_fn_interacted_allparams(theta_c, theta_g, lambda_true);
        [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetac_thetag( moment_fn_interacted, theta_c_grid, theta_g_grid,Z_draws_interacted, alpha, beta, Sigma_conditional);
    else
         moment_fn_interacted = @(theta_c, theta_g) -moment_fn_interacted_allparams(theta_c, theta_g, lambda_true);
        [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetac_thetag( moment_fn_interacted, theta_c_grid, theta_g_grid,Z_draws_interacted, alpha, beta);
    end
    
    interacted_rejection_grids_cell(ds,:) = {grid_lf, grid_rsw, grid_conditional,grid_hybrid};

    %     
%     ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, 'Interacted_Moments/grid', num2str(ds));
%     ds_name = ds_name{:};
%     save( ds_name, 'grid_lf', 'grid_rsw', 'grid_conditional', 'grid_hybrid');
    end

    %Save the cell for the normal moments

    ds_name = strcat( data_output_dir, dirname, 'Basic_Moments/grid_cell');
    ds_name = ds_name{:};
    save( ds_name, 'rejection_grids_cell');
    
    %Save the cell for the interacted moments
    ds_name = strcat( data_output_dir, dirname, 'Interacted_Moments/grid_cell');
    ds_name = ds_name{:};
    save( ds_name, 'interacted_rejection_grids_cell');


end

toc;

%% Aggregate all the results to get the probabilities
for moment_type = {'Basic_Moments/', 'Interacted_Moments/'};
%for moment_type = {'Basic_Moments/'};
    
for dirname = dirnames
    
rejection_prob_lf = zeros( length(theta_c_grid), length(theta_g_grid));
rejection_prob_rsw = rejection_prob_lf;
rejection_prob_conditional = rejection_prob_lf;
rejection_prob_hybrid = rejection_prob_lf;


ds_name = strcat( data_output_dir, dirname, moment_type, 'grid_cell');
ds_name = ds_name{:};
load(ds_name);

for ds = 1:numdatasets

%     ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, moment_type, 'grid', num2str(ds));
%     ds_name = ds_name{:};
%     load(ds_name);
    if( strcmp(moment_type ,'Basic_Moments/'))
    grid_lf = rejection_grids_cell{ds,1};
    grid_rsw = rejection_grids_cell{ds,2};
    grid_conditional = rejection_grids_cell{ds,3};
    grid_hybrid = rejection_grids_cell{ds,4};
    
    else
        grid_lf = interacted_rejection_grids_cell{ds,1};
        grid_rsw = interacted_rejection_grids_cell{ds,2};
        grid_conditional = interacted_rejection_grids_cell{ds,3};
        grid_hybrid = interacted_rejection_grids_cell{ds,4};
    end
    
    rejection_prob_lf = rejection_prob_lf + grid_lf / numdatasets;
    rejection_prob_rsw = rejection_prob_rsw + grid_rsw / numdatasets;
    rejection_prob_conditional = rejection_prob_conditional + grid_conditional / numdatasets;
    rejection_prob_hybrid = rejection_prob_hybrid + grid_hybrid / numdatasets;
    
end
[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);

figure_dir = strcat( figures_output_dir, dirname, moment_type);
figure_dir = figure_dir{:};


contour( xgrid, ygrid, rejection_prob_lf', [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
title('Rejection Probabilities - Least Favorable');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_dir, 'rejection_probs_lf'), 'epsc');

contour( xgrid, ygrid, rejection_prob_rsw', [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
title('Rejection Probabilities - RSW');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_dir, 'rejection_probs_rsw'), 'epsc');

contour( xgrid, ygrid, rejection_prob_hybrid', [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
title('Rejection Probabilities - Hybrid');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_dir, 'rejection_probs_hybrid'), 'epsc');

contour( xgrid, ygrid, rejection_prob_conditional', [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');
title('Rejection Probabilities - Conditional');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_dir, 'rejection_probs_conditional'), 'epsc');


% contourf( xgrid, ygrid, rejection_prob_conditional);
% contour( xgrid, ygrid, rejection_prob_conditional, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on');

end
end
toc;


%% Do a grid just for theta_g using the plus/minus weight moments

% theta_g_grid = -450:50:450;
% 
% numdatasets = 500;
% 
% Z_draws_basic_thetag = Z_draws_basic(5:6,:);
% Z_draws_interacted_thetag = Z_draws_interacted([5:6,11:12,17:18] ,:);
% 
% %Create an inline function to subset columns
% subset_cols = @(mat, cols) mat(:, cols);
% 
% for dirname = dirnames
% 
%     for ds = 1:numdatasets
%     
%     input_dir = strcat( '../../Output/Simulated_Data/', dirname)
%     ds
%     
%     [moment_fn_allparams,moment_fn_interacted_allparams] = generate_moment_fn( input_dir, ds); 
%     
%     %The basic moment fn is columns 5 and 6 of the basic moment fn of all
%     %the parameters (since this doesn't depend on theta_g, i set it to 0)
%     moment_fn = @(theta_g) - subset_cols( moment_fn_allparams(0, theta_g, lambda_true), 5:6);
%     
%     %Compute and save the grids for the basic moments
%     [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetag_only( moment_fn, theta_g_grid,Z_draws_basic_thetag, alpha, beta);
%  
%     ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, 'Basic_Moments/theta_g_grid', num2str(ds));
%     ds_name = ds_name{:};
%     save( ds_name, 'grid_lf', 'grid_rsw', 'grid_conditional', 'grid_hybrid');
%     
%     %Compute and save the grid for the interacted moments
%     moment_fn_interacted = @(theta_g) - subset_cols( moment_fn_interacted_allparams(0, theta_g, lambda_true), [5:6,11:12,17:18]);
% 
%     [grid_lf, grid_rsw, grid_conditional,grid_hybrid] = grids_thetag_only( moment_fn_interacted, theta_g_grid,Z_draws_interacted_thetag, alpha, beta);
%     
%     ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, 'Interacted_Moments/theta_g_grid', num2str(ds));
%     ds_name = ds_name{:};
%     save( ds_name, 'grid_lf', 'grid_rsw', 'grid_conditional', 'grid_hybrid');
%     end
% 
% end
% 
% 
% %% Aggregate the results to get the probabilities for theta_g
% for moment_type = {'Basic_Moments/', 'Interacted_Moments/'};
%     
% for dirname = dirnames
%     
% rejection_prob_lf = zeros( size(theta_g_grid) );
% rejection_prob_rsw = rejection_prob_lf;
% rejection_prob_conditional = rejection_prob_lf;
% rejection_prob_hybrid = rejection_prob_lf;
% 
% for ds = 1:numdatasets
% 
%     ds_name = strcat( '../../Output/Rejection_Grids/Lambda_Constant/', dirname, moment_type, 'theta_g_grid', num2str(ds));
%     ds_name = ds_name{:};
%     load(ds_name);
%     
%     rejection_prob_lf = rejection_prob_lf + grid_lf / numdatasets;
%     rejection_prob_rsw = rejection_prob_rsw + grid_rsw / numdatasets;
%     rejection_prob_conditional = rejection_prob_conditional + grid_conditional / numdatasets;
%     rejection_prob_hybrid = rejection_prob_hybrid + grid_hybrid / numdatasets;
%     
% end
% 
% 
% figure_dir = strcat( '../../Output/Figures/Rejection_Grids/Theta_g_Only/', dirname, moment_type);
% figure_dir = figure_dir{:};
% mkdir(figure_dir);
% 
% plot( repmat( theta_g_grid', 1, 4), [rejection_prob_lf', rejection_prob_rsw', rejection_prob_conditional', rejection_prob_hybrid'] );
% legend( 'LF', 'RSW', 'Conditional', 'Hybrid', 'Location','eastoutside' );
% ylabel( 'Rejection Probability');
% xlabel( 'theta\_g');
% ylim([0, 1]);
% set(gca, 'xtick', [-500:100:500]);
% 
% saveas( gcf, strcat(figure_dir, 'rejection_probs'), 'epsc');
% 
% end
% end
% 
