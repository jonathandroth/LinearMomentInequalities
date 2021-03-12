function [AS_CI, KMS_output] = projected_AS_or_KMS(y_T, X_T, T, Sigma, l, W, method) 

% This function computes projected AS intervals for moment inequalities of
% the form E[Y_t - X_t delta| X_t ] >=0, where the goal is inference on l'delta

% The inputs are:
% y_T = 1/sqrt(T) * sample mean of Y_t
% X_T = 1/sqrt(T) * sample mean of X_t
% T = number of observations
% Sigma = estimate of the cond'l variance E[ Var(Y_t | X_) ]
% l = a vector so that the target parameter is l'delta


%% Set up KMS options

%method      = 'AS';    % Method - either AS or KMS
%method      = 'KMS';    % Method - either AS or KMS
alpha       = 0.05;     % Significance level
KMSopts  = KMSoptions();

KMSopts.CVXGEN = 1;
%KMSopts.CVXGEN_name = [];
KMSopts.DGP = -1; %this tells the functions moments_w, moments_theta, moments_stdev, and moments_gradient what dgp we're using
KMSopts.n = T;
%CVXGEN_name = [];
%CVXGEN_name = "six_moments_two_parameters";

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

%bounds on the parameter space are set to [-10^5,10^5]
LB_theta = repmat(-10^5,size(X_T,2),1); 
UB_theta = repmat(10^5,size(X_T,2),1);

%Store X in KMSoption
KMSopts.X = 1/sqrt(T)* X_T;

%Turn off parallelization
KMSopts.parallel = 0;

KMSopts.S = 0;

num_moments = size(X_T,1);
num_params = size(X_T,2);

if(num_moments == 6 && num_params == 2)
    CVXGEN_name = 'two_parameters_6_moments';
elseif(num_moments == 14 && num_params == 2)
    CVXGEN_name = 'two_parameters_14_moments';
elseif(num_moments == 14 && num_params == 4)
    CVXGEN_name = 'four_parameters_14_moments';
elseif(num_moments == 38 && num_params == 4)
    CVXGEN_name = 'four_parameters_38_moments';
elseif(num_moments == 38 && num_params == 10)
    CVXGEN_name = 'ten_parameters_38_moments';
elseif(num_moments == 110 && num_params == 10)
    CVXGEN_name = 'ten_parameters_110_moments';
else
    warning('Can only use CVXGEN if num_moments and num_parameters match a one of our simulation specs. Using CVX');
    KMSopts.CVXGEN = 0;
    CVXGEN_name = [];
end        

%Increase the # of iterations for KMS
KMSopts.EAM_maxit=1000;

%Compute the AS confidence interval
[AS_CI,KMS_output] = KMS_0_Main(W,theta_0,...
    l,[],LB_theta,UB_theta,[],[],alpha,type,method,kappa,phi,CVXGEN_name,KMSopts);


end
