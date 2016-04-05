ds =1;
load( strcat('../../Output/Simulated_Data/ds', num2str(ds), '.mat' ));

C_moment1 = (J_t_vec == 1) & (J_tminus1_vec == 1);

sum( Pi_vec .* C_moment1) ./ sum( C_moment1) - lambda * ( theta_c + theta_g * mean(g_vec) )