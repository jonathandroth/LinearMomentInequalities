%This function provides the value of eta that solves

% min_{eta,delta} eta

%s.t. Y_T - X_T * delta <= eta

% for the provided y_T and x_T

% This can be compared against a critical value to determine whether there
% we can reject the null hypothesis that E[Y_T - X_T delta] <=0 for some
% delta

function eta_star = test_delta_lp_fn( y_T, X_T)


%Run the linear program

% max_{eta , delta} eta 
% s.t. Y_T -X_T delta <= repmat(eta

%To do this, need to put into the form

% max f X
% st A X <= B

k_delta = size(X_T,2);

f = [1; zeros(k_delta, 1) ];

A_eta = [ ones( size(X_T,1) , 1), zeros(size(X_T,1) , k_delta) ];
A_delta = X_T * [zeros(k_delta,1), eye(k_delta)] ;


A = - (A_eta + A_delta);

b = - y_T;

[~,eta_star] = linprog( f, A, b);
%eta_star = eta_star(1);

end


