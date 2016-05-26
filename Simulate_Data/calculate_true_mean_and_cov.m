%This script calculates the mean and covariance of the moments across the
%500 simulated datasets at the true parameters

%dirnames = { 'Calibrated_SigmaZeta/', 'Calibrated_SigmaZeta_Over4/', 'SigmaZeta_Equal0/'};
%dirnames = { 'Calibrated_SigmaZeta/'};
dirnames = { 'SigmaZeta_Equal0/'};
 
theta_c_true = 129.73;
theta_g_true = -21.38;
lambda_true = 0.386;

g_vec = linspace(12.7, 54.277, 22)'/10; %Divide by 10 to put into 10k lb units
delta_g = g_vec(2) - g_vec(1);

numdatasets = 500;
moments_mat = NaN(numdatasets, 6);

for dirname = dirnames

    for ds = 1:numdatasets
    
    input_dir = strcat( '../../Output/Simulated_Data/', dirname);
    
    ds
    
    [moment_fn_allparams,moment_fn_interacted_allparams] = generate_moment_fn( input_dir, ds); 
    
    %Compute and save the grids for the basic moments
    moment_fn = @(theta_c, theta_g) -moment_fn_allparams(theta_c, theta_g, lambda_true);
    
    moment_ds = mean( moment_fn(theta_c_true, theta_g_true),1 );
    
    moments_mat(ds,:) = moment_ds;
    end
end

mean_moments = mean( moments_mat)
cov_moments = cov( moments_mat)

sd_upper = sqrt( cov_moments(5,5) );
mean_upper = mean_moments(5);

mean_deltapi_upper = - mean_upper + theta_g_true * delta_g;

k = 2;
(mean_deltapi_upper + k * sd_upper) / delta_g


