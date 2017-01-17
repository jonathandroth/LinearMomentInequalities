
lambda_true = .386

lp_create_moments_and_covariances

X_T = X_T_cell{1};
y_T = y_T_cell{1};
Sigma = Sigma_conditional_cell{1};


lp_conditional_test_fn( y_T, X_T, Sigma, alpha)

%Add in 1/100 times the average moment, so that solution to LP is unique
c = 0.01;
k = size(X_T,1);
A = ( eye(k,k) + c*ones(k,k) );

X_T = A * X_T;
y_T = A * y_T;
Sigma = A * Sigma * A';


lp_conditional_test_fn( y_T, X_T, Sigma, alpha)