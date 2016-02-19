%% Test the correlations between the pi_stars

%preferably set J and T to be large before calculating Pi_star_array

pi_star2 = permute( Pi_star_array, [3 1 2] );
cov_mat = zeros(F);

for f = 1:J
cov_mat = cov_mat + cov(pi_star2(:,:,f));
end

cov_mat = 1/J * cov_mat




%% Check convergence

T = 1000;

sigma_nu = 100;
sigma_epsilon = 25;


J_t_startat1 = calculate_offerings( sigma_nu, ...
    sigma_epsilon, ones(F,J) , Epsilon_shocks_array, Eta_shocks_array, Zetaj_shocks_array,Zetajft_shocks_array, Mu_f_array, G_array, F_array );

J_t_startat0 = calculate_offerings( sigma_nu, ...
    sigma_epsilon, zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, Zetaj_shocks_array,Zetajft_shocks_array, Mu_f_array, G_array, F_array );


mean(J_t_startat1(:))
mean(J_t_startat0(:))


meanyrt = @(J,t) mean( mean(J(:,:,t)) );

meanyr_J0 = @(t) meanyrt(J_t_startat0, t);
meanyr_J1 = @(t) meanyrt(J_t_startat1, t);

J0_means = arrayfun( meanyr_J0 , (1:T)' );
J1_means = arrayfun( meanyr_J1 , (1:T)' );

plot( [(1:T)', (1:T)'], [J0_means, J1_means] )


