clear
% This script attempts to replicate Isaiah's power calcs for the basic
% methods, as a check that these are working correctly

alpha = 0.05;
beta = 0.005;

k = 2;
Sigma = eye(k);
%mu_grid = linspace( -10, 10,201);
mu_grid = 0;

grid_length = size(mu_grid, 2);

rejection_prob_lf = zeros( grid_length);
rejection_prob_rsw = rejection_prob_lf;

numsims = 1000;

%Create a matrix of normal draws with mean 0 and covariace sigma
%Each column is a draw
rng(0);
norm_draws = mvnrnd( zeros(1,k) , Sigma, numsims)';

%Create a matrix of standard normals of size k x 10000 for simulating
%critical values
Z_draws = randn(k, 10000);

for s = 1:numsims
    
    z_draw = norm_draws(:, s); 
    
for mu1_index = 1:grid_length
    for mu2_index = 1:grid_length
        
        mu1 = mu_grid(mu1_index);
        mu2 = mu_grid(mu2_index);
        
        %Construct g_T using the norm_draw and the given mus
        g_T = z_draw + [mu1 ; mu2 ];
        
        %Calculate the normalized maximum (used for multiple test)
        R_gt = R_gt_sigma( g_T, Sigma);
        
        %Do the LF test
        cutoff_lf = c_lf(Sigma, alpha, Z_draws);
        test_lf =  ( R_gt > cutoff_lf );
        
        %Do the RSW test
        [cutoff_RSW , mu_tilde] = c_RSW( g_T, Sigma, alpha, beta, Z_draws);
        test_rsw = ( R_gt > cutoff_RSW ) .* max( mu_tilde == 0) ;

        
        rejection_prob_lf(mu1_index, mu2_index) =  rejection_prob_lf(mu1_index, mu2_index) + test_lf;
        rejection_prob_rsw(mu1_index, mu2_index) =  rejection_prob_rsw(mu1_index, mu2_index) + test_rsw;
    end
    
end

end

rejection_prob_lf = rejection_prob_lf / numsims;
rejection_prob_rsw = rejection_prob_rsw / numsims;

rejection_prob_lf
rejection_prob_rsw