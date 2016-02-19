%% Test the correlations between the pi_stars

%preferably set J and T to be large before calculating Pi_star_array

pi_star2 = permute( Pi_star_array, [3 1 2] );
cov_mat = zeros(F);

for f = 1:J
cov_mat = cov_mat + cov(pi_star2(:,:,f));
end

cov_mat = 1/J * cov_mat




