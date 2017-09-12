function [vlo,vup,eta, gamma_tilde] = lp_dual_fn( y_T, X_T,eta, gamma_tilde, Sigma)

%% Calculate vlo and vup
W_T = [ ones(size(X_T,1) , 1) , X_T];
s_T =  ( eye( size(y_T,1) ) - (Sigma * (gamma_tilde * gamma_tilde')) ./ (gamma_tilde' * Sigma * gamma_tilde) ) * y_T;
[vlo, vup] = vlo_vup_dual_fn(eta,s_T,gamma_tilde, Sigma, W_T);

end