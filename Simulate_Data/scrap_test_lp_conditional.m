%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/Single_Thetac/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/Single_Thetac/';


%Specify the groups of firms that have different coefficients
F_group_cell = {[1;2;3;4;5;6;7;8;9]};


%num_F_groups = size(F_group_cell,1);
%l = [1; zeros(num_F_groups-1,1); mean_g];

%Run the main script (this does confidence sets for the mean weight)
xlim_graph = [-10;150];

use_basic_moments = 0;
%% Specify parameters
lp_set_parameters
%% Import data and calculate moments and covariance matrices
lp_create_moments_and_covariances

numdatasets = 3;


lambda_vec = (0.01:1:5.01)';
test1 = NaN(numdatasets, size(lambda_vec,1));
test2 = NaN(size(test1));
i = 1;


for lambda = lambda_vec'
%lambda_true = .386

lambda
lp_create_moments_and_covariances

for ds = 1:numdatasets
X_T = X_T_cell{ds};
y_T = y_T_cell{ds};
Sigma = Sigma_conditional_cell{ds};


test1(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);

%Add in 1/1000 times the average moment, so that solution to LP is unique
c = 0.001;
k = size(X_T,1);
A = ( eye(k,k) + c*ones(k,k) );

X_T = A * X_T;
y_T = A * y_T;
Sigma = A * Sigma * A';


test2(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);


end

i = i+1;
end

test2