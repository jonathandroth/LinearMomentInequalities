clear
% This script attempts to replicate Isaiah's power calcs for the basic
% methods, as a check that these are working correctly

alpha = 0.05;
beta = 0.005;

k = 2;
Sigma = eye(k);
mu_grid = -10:.5:10;
%mu_grid = 0;

grid_length = size(mu_grid, 2);

rejection_prob_lf = zeros( grid_length);
rejection_prob_rsw = rejection_prob_lf;
rejection_prob_conditional = rejection_prob_lf;

%numsims = 1000;
numsims = 100;

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
        
        %Do the conditional test
        %cutoff_conditional = c_conditional(g_T, Sigma, alpha);
        %test_conditional = (R_gt > cutoff_conditional);
        
        T_conditional = c_conditional(g_T, Sigma, alpha);
        test_conditional = T_conditional < alpha;
        
        rejection_prob_lf(mu1_index, mu2_index) =  rejection_prob_lf(mu1_index, mu2_index) + test_lf;
        rejection_prob_rsw(mu1_index, mu2_index) =  rejection_prob_rsw(mu1_index, mu2_index) + test_rsw;
        rejection_prob_conditional(mu1_index, mu2_index) =  rejection_prob_conditional(mu1_index, mu2_index) + test_conditional;
    end
    
end

end

rejection_prob_lf = rejection_prob_lf / numsims;
rejection_prob_rsw = rejection_prob_rsw / numsims;
rejection_prob_conditional = rejection_prob_conditional / numsims;

%Make the contour plots
[ xgrid, ygrid] = meshgrid( mu_grid , mu_grid);
contour(xgrid, ygrid, rejection_prob_lf, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on')
saveas( gcf, '../../Output/power_calcs_lf', 'epsc');

contour(xgrid, ygrid, rejection_prob_rsw, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on')
saveas( gcf, '../../Output/power_calcs_rsw', 'epsc');

contour(xgrid, ygrid, rejection_prob_conditional, [0.05, 0.2, 0.5,0.7,0.9,0.99], 'ShowText', 'on')
saveas( gcf, '../../Output/power_calcs_conditional', 'epsc');
