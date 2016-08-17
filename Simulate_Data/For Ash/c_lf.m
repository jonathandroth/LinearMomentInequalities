
%This function takes as input:
% Sigma: a covariance matrix
% alpha: the significance level
% Z_draws: a matrix of normal draws used to calculate the critical value

%It returns:
% c_lf_alpha: the least favorable critical value with significance alpha
%(assuming a mu of all zeros)

% Zeta_b_mat: the simulated zeta_b's (using Isaiah's notation). These will
% be used in the estimation for the other critical values, such as RSW.

function [c_lf_alpha , Zeta_b_mat ] = c_lf(Sigma, alpha, Z_draws)


D_sigma_minushalf = diag( sqrt( diag(Sigma) ).^(-1) );
Sigma_sqrt = Sigma^(1/2);

%Compute the LF critical values

Zeta_b_mat = D_sigma_minushalf * Sigma_sqrt * Z_draws;
R_b = max(Zeta_b_mat ,[], 1); %Columnwise max of Zeta_b_mat
c_lf_alpha = quantile( R_b, 1-alpha); %take the 1-alpha quantile of R_b


end

