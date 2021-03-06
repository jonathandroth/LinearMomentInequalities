
function [eta, delta, gamma, error_flag] = etahat_fn( y, X, Sigma, varargin)
%% Compute the test statistic etahat used in Andrews, Roth and Pakes for testing moments of the form E[Y_i - X_i delta | Z_i] <= 0
%Inputs
% y: the (scaled) sample average of Y_i (a k x 1 vector)
% X: the (scaled) sample average of X_i (a k x m vector)
% Sigma: the (scaled) estimate of E[ Var(Y_i|Z_i) ] 
%Outputs:
% eta: thet test statistic
% delta: the optimal value of delta in the linear program to compute eta
% gamma: the lagrange multipliers in the program to compute eta
% error_flag: equal to 1 if there was a problem with the LP to compute eta
%Notes
% If y and X are sample averages, Sigma should be E[ Var(Y|X) ] /N
% If Sigma is E[ Var(Y|X) ], then y and X should be averages scaled by sqrt(N)
%The linear program for the test stat is
% min_{eta , delta} eta 
% s.t. Y_T -X_T delta <= eta * sqrt(diag(Sigma)) 



%To do this, need to put into the form

% max f X
% st A X <= B

error_flag =0;

sigma = sqrt(diag(Sigma));

k_delta = size(X,2);

f = [1; zeros(k_delta, 1) ];


A = - [sigma, X]; %This is equivalent to the three lines above


b = - y;

%If options provided, use those; otherwise, set default
if( isempty(varargin) == 0)
    options = varargin{1};
else
    options = optimoptions('linprog','Algorithm','interior-point', 'MaxIter', 20000,'Display', 'off');
end

[delta,eta,flag,~,gamma] = linprog( f, A, b ,[],[],[],[],[], options);

delta = delta(2:end); %Remove eta, which is the first element


if(isstruct(gamma) )
    gamma = gamma.ineqlin;% we put this in an if statement since lambda not defined if there are errors in lp 
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


