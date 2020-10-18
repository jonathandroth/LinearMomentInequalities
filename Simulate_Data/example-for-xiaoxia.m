
alpha = 0.05;
%Need to draw the Z_draws here
Z_draws_interacted = randn(18, 1000)

%NORMALIZE ALL MOMENTS AND Sigma so that all moments have variance 1
Sigma = eye(length(y_T))

%Calculate test statistic
eta = test_delta_lp_fn( y_T, X_T, Sigma);

%Calculate critical value for LFP test
c_alpha = c_lf(Sigma, alpha, Z_draws_interacted);

%Calculate critical value for the 
[c_lp_alpha, eta_draws] = c_lf_lp(X_T,Z_draws_interacted,Sigma,alpha);


lf_test_original(ds,i) = eta > c_alpha;
lf_test_modified(ds,i) = eta > c_lp_alpha;

conditional_test(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);
hybrid_test(ds,i) = lp_hybrid_test_fn( y_T, X_T, Sigma, alpha, alpha/10, eta_draws);