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

function [eta, delta, lambda, error_flag] = etahat_fn( y_T, X_T, Sigma, varargin)


%Run the linear program

% min_{eta , delta} eta 
% s.t. Y_T -X_T delta <= eta * sqrt(diag(Sigma)) 

%To do this, need to put into the form

% max f X
% st A X <= B

error_flag =0;

sigma = sqrt(diag(Sigma));

k_delta = size(X_T,2);

f = [1; zeros(k_delta, 1) ];


A = - [sigma, X_T]; %This is equivalent to the three lines above


b = - y_T;

%If options provided, use those; otherwise, set default
if( isempty(varargin) == 0)
    options = varargin{1};
else
    options = optimoptions('linprog','Algorithm','interior-point', 'MaxIter', 20000,'Display', 'off');
end

[delta,eta,flag,~,lambda] = linprog( f, A, b ,[],[],[],[],[], options);

delta = delta(2:end); %Remove eta, which is the first element


if(isstruct(lambda) )
    lambda = lambda.ineqlin;% we put this in an if statement since lambda not defined if there are errors in lp 
end


%Deal with infinite case
if(flag == -3)
    eta = Inf;
    error_flag = 2;
    warning('Warning: linear program diverged. Setting eta to inf');
    return;
%Deal with other errors in LP
elseif( flag <= 0)
    eta = Inf;
    warning(strcat('Warning: error linear program (flag ', num2str(flag),'). Setting eta to inf'));
    error_flag = 1;
    return;
end



end

