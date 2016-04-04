
%This function returns the list favorable cutoff with signaficance alpha
%and covariance Sigma. It calculates this value using the provided standard normal
%draws (z_draws)

function [c_lf_alpha ] = c_lf(Sigma, alpha, Z_draws)


one_over_sigma_ii = diag( sqrt( diag(Sigma) ) );
Sigma_sqrt = Sigma^(1/2);

%Compute the LF critical values
Zeta_b_mat = one_over_sigma_ii * Sigma_sqrt * Z_draws;
R_b = max(Zeta_b_mat ,[], 1); 
c_lf_alpha = quantile( R_b, 1-alpha);


end

