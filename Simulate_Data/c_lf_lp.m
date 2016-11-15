function [c_lf_lp, eta_vec] = c_lf_lp(X_T, Z_draws, Sigma, alpha)


D_sigma_minushalf = diag( sqrt( diag(Sigma) ).^(-1) );
Sigma_sqrt = Sigma^(1/2);
eps_draws = D_sigma_minushalf * Sigma_sqrt * Z_draws;


numsims = size(Z_draws,2);
eta_vec = NaN(numsims,1);
for s = 1:numsims
    
    y_T_s = eps_draws(:,s);
    eta_vec(s,1) = test_delta_lp_fn(y_T_s, X_T);
end


c_lf_lp = quantile( eta_vec, 1-alpha);

end