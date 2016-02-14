
function sample_moments = sample_moments(Delta_pi_mat,t, g, lambda, theta_c, theta_g) 

%Right now, let's sum across years, firms and products, because dont have
%moments for each product in each year

%Average across products and firms in a given year; average across years


Moments_mat = NaN(size(Delta_pi_mat) );

Moments_mat(:,1) = Delta_pi_mat(:,1) - lambda * theta_c - lambda * theta_g * g;
Moments_mat(:,2) = Delta_pi_mat(:,2) + lambda * theta_c +  lambda * theta_g * g;
Moments_mat(:,3) = Delta_pi_mat(:,3) - theta_c -  theta_g * g;
Moments_mat(:,4) = Delta_pi_mat(:,4) + theta_c +  theta_g * g;
Moments_mat(:,5) = Delta_pi_mat(:,5) - theta_g;
Moments_mat(:,6) = Delta_pi_mat(:,6) + theta_g;

sample_moments = mean( grpstats( Moments_mat, t, {'nanmean'} ) );

end