lambda_vec = 0.01:.5:5.01
test1 = NaN(size(lambda_vec));
test2 = NaN(size(lambda_vec));
i = 1;

for lambda_true = lambda_vec
%lambda_true = .386

lambda_true
lp_create_moments_and_covariances

X_T = X_T_cell{1};
y_T = y_T_cell{1};
Sigma = Sigma_conditional_cell{1};


test1(i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);

%Add in 1/1000 times the average moment, so that solution to LP is unique
c = 0.001;
k = size(X_T,1);
A = ( eye(k,k) + c*ones(k,k) );

X_T = A * X_T;
y_T = A * y_T;
Sigma = A * Sigma * A';


test2(i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);

i = i+1;
end

[test1', test2']