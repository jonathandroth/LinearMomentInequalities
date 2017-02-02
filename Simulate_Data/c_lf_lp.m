%This function computes the (modified) least favorable critical values for
%the linear programming approach

%It takes as input

%X_T: the matrix of derivatives, i.e. true moments are of hte form Y_T -
%  X_T * delta

%Z_draws: a matrix of standard normal draws

%Sigma: a covariance matrix

%alpha: significance level

%Note Sigma and X_T are assumed to have been normalized so that Sigma has
%diagonal of ones



function [c_lf_lp, eta_vec] = c_lf_lp(X_T, Z_draws, Sigma, alpha)

%Note: x_T and Sigma are assumed to be rescaled so that Sigma has diagnol
%of 1
%Rescale the covariance matrix to have diagnol of 1
    % X_T is assumed to have already been re-scaled
    %Note: this is sort of an odd placement, where we do this rescaling
    %inside and the other rescaling outside
%D_sigma_minushalf = diag( sqrt( diag(Sigma) ).^(-1) );
Sigma_sqrt = Sigma^(1/2);
%eps_draws = D_sigma_minushalf * Sigma_sqrt * Z_draws;
eps_draws = Sigma_sqrt * Z_draws;

numsims = size(Z_draws,2);
eta_vec = NaN(numsims,1);
for s = 1:numsims
    
    y_T_s = eps_draws(:,s);
    eta_vec(s,1) = test_delta_lp_fn(y_T_s, X_T);
end


c_lf_lp = quantile( eta_vec, 1-alpha);

end