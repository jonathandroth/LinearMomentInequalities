function [J_t_array, J_tminus1_array, Pi_array] = calculate_offerings( sigma_nu, ...
    sigma_epsilon, J_0 , Epsilon_shocks_array, Eta_shocks_array, Zetaj_shocks_array,Zetajft_shocks_array, Mu_f_array, G_array, F_array )

global T;
global J;
global F;

global lambda;
global theta_c;
global theta_g;



%For now ,we assume that sigma_zetaj = sigma_nu, sigma_zetajft =
%sigma_epsilon
sigma_zetaj = sigma_nu;
sigma_zetajft = sigma_epsilon;

%Calculate the Nu draws and Epsilon draws using the shocks and then given
%means and covariances

%Eta_draws_array =  Mu_f_array + gamma* theta_g * G_array + sigma_nu * Eta_shocks_array;
%Epsilon_draws_array = sigma_epsilon * Epsilon_shocks_array;
%Pi_star_array = Eta_draws_array + Epsilon_draws_array;


J_tminus1_array = NaN(F,J,T);
J_tminus1_array(:,:,1) = J_0;

J_t_array = NaN(F,J,T);
Pi_star_array = NaN(F,J,T);

for t=1:T
   
    J_tminus1 = J_tminus1_array(:,:,t);
    
    
    
    %Compute the eta and epsilon shocks for period t by shifting and
    %scaling the standard normal shocks by the appropriate means and
    %variances
    %Right now, eta has mean mu_f + theta_g g if j is not in J_tminus1; eta
    %has mean mu_f +(lambda * theta_g g) if j is in J_tminus1
    Eta_t = Mu_f_array(:,:,t) + ( J_tminus1 * lambda + (1 - J_tminus1) ) .* (  theta_g * G_array(:, :, t) ) + sigma_nu * Eta_shocks_array(:,:,t) ;
    Epsilon_t = sigma_epsilon * Epsilon_shocks_array(:,:,t);
    
    Pi_star_t = Eta_t + Epsilon_t;
    Pi_star_array(:,:,t) = Pi_star_t;
    
    %Calculate profit (including FCs) from adding J in period t, using fixed costs
    %dependent on J_tminus1
    Pi_t = Pi_star_t - (J_tminus1 .* lambda + (1 - J_tminus1) ) .* (theta_c + theta_g * G_array(:, :, t) );
    
    %j is in J_t if profit including fixedcosts is greater than 0
    J_t = (Pi_t >= 0);
    J_t_array(:,:,t) = J_t;
    
    %Update the J_tminus1 array
    if t< T
    J_tminus1_array(:,:, t+1) = J_t;
    end
end


%Add shocks to pi_star to generate pi
Pi_array = Pi_star_array .* (1 + sigma_zetaj * Zetaj_shocks_array + sigma_zetajft * Zetajft_shocks_array);



end