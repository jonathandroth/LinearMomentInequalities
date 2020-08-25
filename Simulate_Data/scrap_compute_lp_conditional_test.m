clear
clc

load('test.mat')

%%
[reject, eta, delta, lambda, pval] = lp_conditional_test_fn(y_T, X_T, sigma, 0.05);

%%
[eta_star, delta_star, lambda, error_flag] = test_delta_lp_fn( y_T, X_T, sigma, optimoptions('linprog','Algorithm','dual-simplex','TolFun', 10^-8, 'Display', 'off', 'MaxIter', 100000))

%%
[vlo,vup,eta, gamma_tilde] = lp_dual_fn( y_T, X_T, eta, lambda, sigma)

%%
sdVec = sqrt(diag(sigma));
W_T = [ sdVec , X_T];
s_T =  ( eye( size(y_T,1) ) - (sigma * (gamma_tilde * gamma_tilde')) ./ (gamma_tilde' * sigma * gamma_tilde) ) * y_T;
[vlo, vup] = vlo_vup_dual_fn(eta,s_T,gamma_tilde, sigma, W_T);


