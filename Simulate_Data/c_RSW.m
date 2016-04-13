


function c_RSW = c_RSW( g_T, Sigma, alpha, beta, Z_draws)

[c_beta_lf, Zeta_b_mat] = c_lf(Sigma, beta, Z_draws);


mu_tilde =  min( g_T + sqrt( diag( Sigma) ) .* c_beta_lf , 0); 

Zeta_b_mat = Zeta_b_mat + mu_tilde; %Add mu_tilde to the zeta_b_mat from the lf approach (which assumes mu = 0)
R_b = max(Zeta_b_mat ,[], 1); %Columnwise max of Zeta_b_mat
c_RSW = quantile( R_b, 1-alpha+beta); %take the 1-alpha quantile of R_b



