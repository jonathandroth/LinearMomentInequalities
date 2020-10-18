%% Run simulation of APR inference for tau2 treating tau1 as a nuisance
% paramter, where there is 1 pre-period

A = [1, 1, 0; 0, -2, 1] ;
A = [A;-A];

%Create X that treats tau1 as the nuisance parameter
tauCoefs = [0; 1 ; 0];
X = A *tauCoefs;

%betahat = [-1;1;3];
betahat = [0;1;4];
%betahat = [-8;1;6];
%betahat = [-20;1;50];

%Sigma = eye(length(betahat));
%Sigma = [1,0.5,0.5; 0.5, 1, 0.5; 0.5, 0.5, 1];
Sigma = [2,0.5,0.5; 0.5, 1, 0.5; 0.5, 0.5, 1];


M = 0;
tau2 = 0;

y = A * (betahat - [0;0;1]*tau2) - M;

SigmaMoments = A * Sigma * A';
sdMatMoments = diag( sqrt( diag( SigmaMoments ) ) );


[eta, tau1, lambda] = test_delta_lp_fn( sdMatMoments^(-1) * y, sdMatMoments^(-1) *  X)


%% Instead remove the nuisance parameter and use the moment |(betahat2 -
%%% tau2) + 2betahat-1| < 3M

A_nonuisance = [2,0,1];
A_nonuisance = [A_nonuisance;-A_nonuisance];

y_nonuisance = A_nonuisance * (betahat - [0;0;1]*tau2) - 3*M;
X_nonuisance = zeros(size(y_nonuisance));


SigmaMoments_nonuisance = A_nonuisance * Sigma * A_nonuisance';
sdMatMoments_nonuisance = diag( sqrt( diag( SigmaMoments_nonuisance ) ) );
[eta_nonuisance, ~, lambda_nonuisance] = test_delta_lp_fn(y_nonuisance, X_nonuisance)


%% Compare results
(lambda' * A) ./ (lambda_nonuisance' * A_nonuisance)

A * (betahat - [0;0;1]*tau2 - [0;1;0]*tau1)
A * (betahat - [0;0;1]*tau2 )

betahat - [0;0;1]*tau2