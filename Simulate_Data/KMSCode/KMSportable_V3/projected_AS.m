function [AS_CI, KMS_output] = projected_AS(y_T, X_T, T, Sigma, l, W) 

% This function computes projected AS intervals for moment inequalities of
% the form E[Y_t - X_t delta| X_t ] >=0, where the goal is inference on l'delta

% The inputs are:
% y_T = 1/sqrt(T) * sample mean of Y_t
% X_T = 1/sqrt(T) * sample mean of X_t
% T = number of observations
% Sigma = estimate of the cond'l variance E[ Var(Y_t | X_) ]
% l = a vector so that the target parameter is l'delta


%% Set up KMS options

method      = 'AS';    % Method - either AS or KMS
alpha       = 0.05;     % Significance level
KMSopts  = KMSoptions();

KMSopts.CVXGEN = 0;
KMSopts.CVXGEN_name = [];
KMSopts.DGP = -1; %this tells the functions moments_w, moments_theta, moments_stdev, and moments_gradient what dgp we're using
KMSopts.n = T;
CVXGEN_name = [];

%% Parameters
type = 'two-sided';       % Two-sided or one sided test?  Set to 'one-sided-UB' or 'one-sided-LB' or 'two-sided'
kappa = NaN;              % Default kappa function
phi   = NaN;              % Default GMS function      


% Draw T normal draws from a Normal with mean yBar and variance Sigma
%This makes the internal bootstrap in KMS code approximate the cond'l
%distribution
if(isnan(W))
yBar = 1/sqrt(T) * y_T;
rng(0); 
W = mvnrnd(yBar, Sigma, T);
end

theta_0 = zeros(size(X_T,2),1);%inital guess is theta =0

%bounds on the parameter space are set to [-100,100]
LB_theta = repmat(-100,size(X_T,2),1); 
UB_theta = repmat(100,size(X_T,2),1);

%Store X in KMSoption
KMSopts.X = X_T;

%Turn off parallelization
KMSopts.parallel = 0;



%Compute the AS confidence interval
[AS_CI,KMS_output] = KMS_0_Main(W,theta_0,...
    l,[],LB_theta,UB_theta,[],[],alpha,type,method,kappa,phi,CVXGEN_name,KMSopts);


end
