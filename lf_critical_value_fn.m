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



function [lf_critical_value_fn, eta_vec] = lf_critical_value_fn(X, Z_draws, Sigma, alpha)


Sigma_sqrt = Sigma^(1/2);
eps_draws = Sigma_sqrt * Z_draws;

numsims = size(Z_draws,2);
eta_vec = NaN(numsims,1);
for s = 1:numsims
    
    y_T_s = eps_draws(:,s);
    eta_vec(s,1) = etahat_fn(y_T_s, X, Sigma);
end


lf_critical_value_fn = quantile( eta_vec, 1-alpha);

end
