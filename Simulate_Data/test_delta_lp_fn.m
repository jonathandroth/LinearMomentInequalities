%This function provides the value of eta that solves

% min_{eta,delta} eta

%s.t. Y_T - X_T * delta <= eta

% for the provided y_T and x_T

% This can be compared against a critical value to determine whether
% we can reject the null hypothesis that E[Y_T - X_T delta] <=0 for some
% delta

%The function takes an optional value to be passed to the linear program

%The function also returns the values
% delta-star = the value of delta at the optimum
% lambda = the lagrange multipliers at the optimum

function [eta_star, delta_star, lambda] = test_delta_lp_fn( y_T, X_T, varargin)


%Run the linear program

% min_{eta , delta} eta 
% s.t. Y_T -X_T delta <= repmat(eta

%To do this, need to put into the form

% max f X
% st A X <= B


k_delta = size(X_T,2);

f = [1; zeros(k_delta, 1) ];

%A_eta = [ ones( size(X_T,1) , 1) , zeros(size(X_T,1) , k_delta) ];
%A_delta = X_T * [zeros(k_delta,1), eye(k_delta)] ;
% A = - (A_eta + A_delta);

A = - [ones( size(X_T,1) , 1), X_T]; %This is equivalent to the three lines above


b = - y_T;

%If options provided, use those; otherwise, set default
if( isempty(varargin) == 0)
    options = varargin{1};
else
    options = optimoptions('linprog','Algorithm','interior-point', 'MaxIter', 2000);
end

[delta_star,eta_star,flag,~,lambda] = linprog( f, A, b ,[],[],[],[],[], options);

if(isstruct(lambda) )
    lambda = lambda.ineqlin;
else
    eta_star = Inf;
    error('Warning: linear program did not converge. Setting eta to inf');
end

delta_star = delta_star(2:end); %Remove eta, which is the first element

%Deal with infinite cases
if(flag == -2 || flag == -4 || flag == -5)
    eta_star = inf;
end

end


