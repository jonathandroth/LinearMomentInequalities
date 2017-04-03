

cov_conditional = conditional_variance_fn(Y_basic, A_basic, diagonal);
cov_unconditional = cov(Y_basic + A_c_basic * theta_c_true + A_g_basic * theta_g_true ) 
%cov_unconditional = cov(Y_basic ) 


cov_conditional
cov_unconditional

[~,p] = chol( cov_unconditional - cov_conditional )

eig(cov_unconditional - cov_conditional )