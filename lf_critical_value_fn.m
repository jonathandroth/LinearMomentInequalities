
function [lf_critical_value, eta_vec] = lf_critical_value_fn(X, Z_draws, Sigma, alpha)
%% Compute the least-favorable critical value from Andrews, Roth and Pakes for testing moments of the form E[Y_i - X_i delta | Z_i] <= 0
%Inputs:
% X: the (scaled) sample average of X_i (a k x m vector)
% Z_draws: a matrix of standard normal draws of dimension k x S, where S is
% the number of simulations for the CV (we recommend 1000). This can be
% drawn once for all tests used for computational stability.
% Sigma: the (scaled) estimate of E[ Var(Y_i|Z_i) ]. (Should be the same as
% Sigma used for calculating the test statistic in the function etahat_fn)
% alpha: the size of the test (e.g. 0.05 for 5% significance)
%Outputs:
% lf_critical_value: the critical value of the test
% eta_vec: the simulated values of eta. The CV is the 1-alpha quantile of
% this vector. This vector can be used as input to the hybrid test for
% computational efficiency

Sigma_sqrt = Sigma^(1/2);
eps_draws = Sigma_sqrt * Z_draws;

numsims = size(Z_draws,2);
eta_vec = NaN(numsims,1);
for s = 1:numsims
    
    y_T_s = eps_draws(:,s);
    eta_vec(s,1) = etahat_fn(y_T_s, X, Sigma);
end


lf_critical_value = quantile( eta_vec, 1-alpha);

end
