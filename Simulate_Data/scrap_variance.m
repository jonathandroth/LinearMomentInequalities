

cov_conditional = conditional_variance_fn(Y_basic, A_basic, diagonal);
%cov_unconditional = cov(Y_basic + A_c_basic * theta_c_true + A_g_basic * theta_g_true ) ;
cov_unconditional = cov(Y_basic + A_c_basic * (theta_c_true) + A_g_basic * (theta_g_true+20) ) ;
%cov_unconditional = cov(Y_basic ) 


cov_conditional
cov_unconditional

[~,p] = chol( cov_unconditional - cov_conditional )
eig(cov_unconditional - cov_conditional )

Z_draws = randn( size(cov_conditional,1) , 10000);

cutoff_conditional = c_lf(cov_conditional, alpha, Z_draws)
cutoff_unconditional = c_lf(cov_unconditional, alpha, Z_draws)

[cutoff_conditional * diag(cov_conditional),...
cutoff_unconditional * diag(cov_unconditional)]

g_T = -sum(Y_basic + A_c_basic * (theta_c_true) + A_g_basic * (theta_g_true+20))' / sqrt( size(Y_basic,1) );

basic_tests(g_T, cov_conditional, Z_draws, alpha, beta)

basic_tests(g_T, cov_unconditional, Z_draws, alpha, beta)
