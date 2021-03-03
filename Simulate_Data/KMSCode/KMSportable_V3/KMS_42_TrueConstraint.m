function [c,ceq] = KMS_42_TrueConstraint(theta,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions)
% This function computes the constraint of the min/max problem on
% Pg 13, Eq 2.16.   We use fmincon to find the upper/lower bounds of the
% confidence set.  This method is likely to perform poorly in realistic
% problems.  This method is used as an error-check in simple problems 
% to ensure that EAM is working correctly.
%
% NB: We do not output gradient of the constraint set, because we do not
% know what the gradient of the critical value C(theta) is.  One method is
% to compute the gradient of the moments analytically, and use numerical
% gradients to approximate dC(theta).

% INPUTS:
%   theta               dim_p-by-1 parameter vector
%
%   f_ineq,f_eq         Empirical moments
%
%   f_ineq_keep,f_eq_keep       Moments to keep
%
%  f_stdev_ineq,f_stdev_eq     Standard deviation of empirical moments
%
%   G_ineq,G_eq         Bootstrapped and recentered moments
%
%   KMSoptions.         This is a structure of additional inputs held
%                       constant over the program.  In the 2x2 entry game,
%                       KMSoptions includes the support for the covariates
%                       and the probability of support point occuring.
%
% OUTPUT:
%   c                   J-by-1 vector of inequality constraints


%% Extract relevant information from KMSoptions
J1          = KMSoptions.J1;
J2          = KMSoptions.J2;
J           = KMSoptions.J;     
n           = KMSoptions.n;

%% Standardized moments
[g_ineq,g_eq] = moments_theta(theta,J1,J2,KMSoptions);

% Standardized momoments
h_theta = sqrt(n)*(([f_ineq;f_eq] + [g_ineq;g_eq])./[f_stdev_ineq;f_stdev_eq]);

% Drop moments close to boundary
f_keep  = [f_ineq_keep;f_eq_keep];
h_theta(f_keep == 0,:) = [];

%% Compute c(theta)
 c_theta = KMS_31_Estep(theta.',f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
 
%% Now compute nonlinear constraint
c = h_theta - c_theta;
ceq = [];

end