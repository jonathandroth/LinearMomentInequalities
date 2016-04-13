
%This function returns the list favorable cutoff with significance alpha
%and covariance Sigma. It calculates this value using the provided standard normal
%draws (Z_draws)

function [c_lf_alpha , Zeta_b_mat ] = c_lf(Sigma, alpha, Z_draws)


D_sigma_minushalf = diag( sqrt( diag(Sigma) ).^(-1) );
Sigma_sqrt = Sigma^(1/2);

%Compute the LF critical values

Zeta_b_mat = D_sigma_minushalf * Sigma_sqrt * Z_draws;
R_b = max(Zeta_b_mat ,[], 1); %Columnwise max of Zeta_b_mat
c_lf_alpha = quantile( R_b, 1-alpha); %take the 1-alpha quantile of R_b


end

