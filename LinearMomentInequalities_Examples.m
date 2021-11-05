% This script provides examples of how to use the functions for the moment
% inequality tests in Andrews, Roth, and Pakes
% We assume that the parameters satisfy E[ Y_i(beta_0) - X_i \delta | Z] >= 0,
% where beta_0 is the target parameter and \delta is a nuisance parameter,
% and X_i is non-random conditional on Z
% The code takes as input y (the average of the Y_i, scaled by sqrt(N)), X (the average of the
% X_i, scaled by sqrt(N)), and Sigma (an estimate of the conditional variance E[ Var(Y_i|Z_i)
% ]
% ]
% Important note: if Sigma is an estimate of E[ Var(Y_i|Z_i)], then y and X
% should be the sample averages *scaled by sqrt(N)*. Alternatively, if y
% and X are sample averages of y_i and X_i, then Sigma should be E[
% Var(Y_i|Z_i) / N. 

% Given a matrix Y where each row corresponds with a realization of Y_i,
% and a matrix Z where each row corresponds with Z_i, one can estimate the
% Sigma using the provided conditional_variance_function, which implements
% the matching method of Abadie Imbens and Zhang (2014) 

% We presume the calculation of Y_i, X_i, and Sigma has been done offline.
% Note that if Y_i or X_i depends on beta0, then Sigma must be calculated
% for each candidate value of beta0 to form a confidence set

% See below for details about the case where beta0 enters linearly, i.e.
% Y(beta0) = Y_i + beta0*X_beta


%Example values of y,X, and Sigma
y = [1;1;0];
X = [1;-1;0];
Sigma = eye(length(y));

%Set significance level
alpha = 0.05; 
%Draw a matrix of standard normals to be used for the LF test
    %It is best to draw this once at the beginning for stability / computational efficiency
    %reasons
Z_draws = randn(1000, length(y))'; %this uses 1000 draws for the LF critical value calc.


%Compute the LF critical value (save the draws, which can be used for the
%hybrid test)
[cv_lf, draws] = lf_critical_value_fn(X,Z_draws, Sigma, alpha);

%Compute the test statistic
etahat = etahat_fn(y,X,Sigma);

%LF test rejects if etahat > cv_lf
reject_lf = eta > cv_lf;

%Run the conditional test
reject_conditional = conditional_test_fn(y,X,Sigma,alpha);

%Run the hybrid test
    %We pass the draws from the LF test to hybrid function. This is
    %optional but can greatly increase speed
reject_hybrid = hybrid_test_fn(y,X,Sigma,alpha, alpha/10, draws)


%If beta_0 enters non-linearly, then the procedure above should be done for
%each candidate values of beta0. Then the values of beta0 that are not
%rejected can be collected to form a confidence set


%If the target parameter enters the moments linearly, there are substantial
%computational savings. First, Sigma need only be calculated once. Second,
%the LF confidence interval can be calculated with a linear program without
%using test inversion. We provide an example below where we are interested
%in the first component of delta

%For the sake of the example, add another column to X so that delta has two dimension.
% Suppose we are interested in the first component of delta
X = [[1;-1;0] , [0;0;0]];

l = [1;0]; %target parameter is l'delta

%Compute the LF critical value
[cv_lf, draws] = lf_critical_value_fn(X,Z_draws, Sigma, alpha);

%Compute the LF confidence interval directly using a linear program
lf_ci_for_linear_params(y,X,Sigma,l,cv_lf) 


%If we want a confidence set for the target parameter for the
%hybrid/conditional test, we need to test inversion. But we get
%computatioanl savings because Sigma need not be recalculated for each
%gridpoint, and we can use the LF results of LF test to reduce computation
%of the hybrid's first stage

beta0_grid = -3:.1:3; %grid for target parameter beta0 = l'delta

conditional_test_grid = NaN(length(beta0_grid),1);
hybrid_test_grid = NaN(length(beta0_grid),1);


%Change of basis so that first component of X corresponds with the
%direction l (not needed if l = (1,0,...0)'
M = [l' ; null( repmat(l',size(l)) )'];
X = X * M^(-1);

Xtilde = X(:,2:end); %Create Xtilde, which corresponds with the component of X orthogonal to l
X_l = X(:,1); %Create x_l which corresponds with the component of X in direction l


for s = 1:length(beta0_grid)
    
   beta0 = beta0_grid(s);
   y_s = y - X_l * beta0;
   
   conditional_test_grid(s) = conditional_test_fn(y_s,Xtilde,Sigma,alpha);
   hybrid_test_grid(s) = hybrid_test_fn(y_s,Xtilde,Sigma,alpha,alpha/10,draws);
end

%Confidence set is points where the test doesn't reject
conditional_CS = beta0_grid( conditional_test_grid == 0) 
hybrid_CS = beta0_grid( hybrid_test_grid == 0) 
