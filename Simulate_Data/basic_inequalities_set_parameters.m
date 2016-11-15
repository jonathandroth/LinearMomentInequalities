
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
nummoments_interacted = 30;

%Create a matrix of standard normals of size k x 10000 for simulating
%critical values
rng(0);
Z_draws_interacted = randn(nummoments_interacted, 10000);
Z_draws_basic = Z_draws_interacted( 1:nummoments_basic,:);

numdatasets = 500;

nummarkets = 500; %This is the number of markets to sample from the long chain

%dirnames = { 'Calibrated_SigmaZeta/'};
 