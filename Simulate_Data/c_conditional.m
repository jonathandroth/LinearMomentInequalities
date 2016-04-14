

function T = c_conditional(g_T, Sigma, alpha)

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

T = 1 -  ( normcdf(h_Ti) - normcdf(v_lo) ) ./ (1 - normcdf(v_lo) );

t_i = h_Ti;
l_i = v_lo ./ sqrt( sigma_diag(is_max == 1) );
u_i = inf;

T_bound = 1- (  (sqrt( 4 + l_i^2) + l_i)^(-1) -  (sqrt( 2 + t_i^2) + t_i)^(-1) * exp(-0.5 *(t_i^2 - l_i^2) ) ) ...
    / ( ( sqrt( 2 + l_i^2) + l_i)^(-1) - (sqrt( 4 + u_i^2) + u_i)^(-1) * exp(-0.5 * (u_i^2 - l_i^2)) );

T_bound( min(T_bound, h_Ti) < 0 ) = 1;


T = min( T, T_bound);

% %Use this to get a cutoff
% zeta_lo = normcdf( v_lo ./ sqrt( sigma_diag( is_max == 1) )  ); 
% c_cond = norminv( 1 - alpha + alpha * zeta_lo);


end