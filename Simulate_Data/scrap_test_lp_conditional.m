%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/Single_Thetac/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/Single_Thetac/';

filename_graph = 'lambda_rejection_probabilities';
%Specify the groups of firms that have different coefficients
F_group_cell = {[1;2;3;4;5;6;7;8;9]};


%num_F_groups = size(F_group_cell,1);
%l = [1; zeros(num_F_groups-1,1); mean_g];

use_basic_moments = 0;
%% Specify parameters
lp_set_parameters
%% Import data and calculate moments and covariance matrices
lp_create_moments_and_covariances


lambda_vec = (0.01:.2:5.01)';
conditional_test_noadjustment = NaN(numdatasets, size(lambda_vec,1));
conditional_test_adjustment = NaN(size(conditional_test_noadjustment));
lf_test_original = NaN(size(conditional_test_noadjustment));
lf_test_modified = NaN(size(conditional_test_noadjustment));

i = 1;


for lambda = lambda_vec'
%lambda_true = .386

lambda
lp_create_moments_and_covariances

parfor ds = 1:numdatasets
X_T = X_T_cell{ds};
y_T = y_T_cell{ds};
Sigma = Sigma_conditional_cell{ds};


conditional_test_noadjustment(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);

%Add in 1/1000 times the average moment, so that solution to LP is unique
c = 0.001;
k = size(X_T,1);
A = ( eye(k,k) + c*ones(k,k) );

X_T = A * X_T;
y_T = A * y_T;
Sigma = A * Sigma * A';


conditional_test_adjustment(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);

%Do non-conditional tests
c_alpha = c_lf(Sigma, alpha, Z_draws_interacted);
c_lp_alpha = c_lf_lp(X_T,Z_draws_interacted(:,1:numsims_lp),Sigma,alpha);
eta = test_delta_lp_fn( y_T, X_T);

lf_test_original(ds,i) = eta > c_alpha;
lf_test_modified(ds,i) = eta > c_lp_alpha;

end

i = i+1;
end

%% Save the results

    ds_dir = strcat( data_output_dir, 'Interacted_Moments/');
    mkdir(ds_dir);
    save( strcat(ds_dir, 'lambda_results'), 'conditional_test_adjustment', 'conditional_test_noadjustment', 'lf_test_original', 'lf_test_modified');


%% Plot
plot(lambda_vec, [mean(conditional_test_adjustment);...
                  mean(conditional_test_noadjustment);...
                  mean(lf_test_original);...
                  mean(lf_test_modified)] ) 
              
              
legend( 'Conditional (Unshifted)', 'Conditional (Shifted)', 'LF','LF (modified)', 'Location','eastoutside' );
ylabel('Rejection Probability');
xlabel('Lambda');
title('Rejection Probabilities for Lambda');
saveas( gcf, strcat(figures_output_dir,filename_graph ), 'epsc');

