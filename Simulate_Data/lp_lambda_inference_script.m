%% Specify parameters
lp_set_parameters

%% Import data and calculate moments and covariance matrices and accept/reject over grid
lp_lambda_moments_and_covariances

    

%% Compute identified set

tic;

display('Began identified set calc');

lp_lambda_identified_set
    
display('Finished identified set calc');

toc;


 %% Plot

lp_lambda_plot