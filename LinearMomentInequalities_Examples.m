y = [1;1;0];
X = [1;-1;0];

Sigma = eye(length(y));

Z_draws = randn(1000, length(y))';
alpha = 0.05; 

[c_lf, draws] = lf_critical_value_fn(X,Z_draws, Sigma, alpha);

%XX add something about studentiziation here

eta = etahat_fn(y,X,Sigma)

conditional_test_fn(y,X,Sigma,alpha)
hybrid_test_fn(y,X,Sigma,alpha, alpha/10, draws)

%%If you're using both LF / Hybrid at the same time, can use draws in
%%hybrid XX 

% If you want a CI for a component of delta, you can use this function XX
cs_linear_delta_lp_fn(y,X,Sigma,1,c_lf) 


%If you want hybrid for a linear parameter, then can pass draws
