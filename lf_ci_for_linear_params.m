%This function provides a confidence set for a linear combination of the
%elements of the delta vectors for moments of the form Y_T - X_T delta >= 0

%It takes as inputs:

%y_T, x_T: data for the moments
% l: a CS is provided for the linear combination l' * delta
% c_alpha: the critical value for the moments


%It returns:
% cs: which is a vector of the upper and lower bounds of the CS: [ub, lb] 

function [cs, slack_lb, slack_ub] = lf_ci_for_linear_params( y, X, Sigma, l, lf_cv) 

%% Compute the least-favorable CI from Andrews, Roth and Pakes for testing moments of the form E[Y_i - X_i delta | Z_i] <= 0, in the special case where the target parameter is l'delta
%Inputs
% y: the (scaled) sample average of Y_i (a k x 1 vector)
% X: the (scaled) sample average of X_i (a k x m vector)
% Sigma: the (scaled) estimate of E[ Var(Y_i|Z_i) ]
% l: The target parameter is l'delta
% alpha: the size of the test (e.g. 0.05 for 5% significance)
% lf_cv: The least favorable CV outputted by lf_critical_value_fn
%Outputs:
% cs: a two-dimensional vector with the lower and upper bounds of the
% confidence set
% slack_lb: the slackness of the moments at the lb bound of the cs
% slack_ub: the slackness of the moments at the ub bound of the cs
%Notes
% If y and X are sample averages, Sigma should be E[ Var(Y|X) ] /N
% If Sigma is E[ Var(Y|X) ], then y and X should be averages scaled by sqrt(N)
% Sigma used here should be the same as the sigma used for lf_cv


sigma = sqrt(diag(Sigma));

A = -X;
b = lf_cv * sigma - y;

[delta_lb,lb] = linprog(l , A, b, [], [], [], [],  optimoptions('linprog','TolFun', 10^(-8), 'Display','off')) ; 
[delta_ub,ub] = linprog(-l, A, b, [], [], [], [],  optimoptions('linprog','TolFun', 10^(-8), 'Display','off'));
ub = -ub; %Need to take negative of the minimum to get the max

cs = [lb; ub];

%Compute slack at upper and lower bound solutions
slack_lb = A * delta_lb - b;
slack_ub = A * delta_ub - b;


end
