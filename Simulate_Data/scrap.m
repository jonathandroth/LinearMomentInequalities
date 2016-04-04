
load('../../Output/Simulated_Data/calibrated_params');

 [J_t_array, J_tminus1_array, Pi_star_array] = calculate_offerings( sigma_nu, ...
        sigma_eps, zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, mu_f, G_array );
    
 
 Jt1 =  J_t_array(:,:,1);
 Jt0 =  J_tminus1_array(:,:,1);
 
 Pistar1 = Pi_star_array(:,:,1);
 
 testJ = (Pi_star_array -  (J_tminus1_array * lambda + (1 - J_tminus1_array))...
     .* (theta_c + G_array *theta_g)      ) >= 0 ;
 
 diffJ = J_t_array - testJ;
 
 sum( abs( diffJ(:) ) )
 