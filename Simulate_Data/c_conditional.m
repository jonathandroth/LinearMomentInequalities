

function c_cond = c_conditional(g_T, Sigma, alpha)

max_g = max(g_T);
is_max = (g_T == max_g);



%Select h_T for the non-max entries, and max-entry
h_T = g_T ./ sqrt( diag(Sigma) );
h_Tj = h_T( is_max == 0);
h_Ti = h_T( is_max == 1);

%Costruct omega_ij
sigma_diag = diag(Sigma); 
omega_ij = Sigma( is_max == 0, is_max == 1) ./  ( sqrt( sigma_diag( is_max == 0) ) .* sqrt( sigma_diag(is_max == 1) ) );

v_lo = max( (h_Tj - omega_ij .* h_Ti) ./ ( 1 - omega_ij) ) ;

zeta_lo = normcdf( v_lo ./ sqrt( sigma_diag( is_max == 1) )  ); 

c_cond = norminv( 1 - alpha + alpha * zeta_lo);


end