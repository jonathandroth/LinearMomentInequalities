%This function provides a confidence set for a linear combination of the
%elements of the delta vectors for moments of the form Y_T - X_T delta >= 0

%It takes as inputs:

%y_T, x_T: data for the moments
% l: a CS is provided for the linear combination l' * delta
% c_alpha: the critical value for the moments


%It returns:
% cs: which is a vector of the upper and lower bounds of the CS: [ub, lb] 

function [cs, slack_lb, slack_ub] = cs_linear_delta_lp_fn( y_T, X_T, l, c_alpha) 


A = -X_T;
b = c_alpha - y_T;

[delta_lb,lb] = linprog(l , A, b, [], [], [], [],  optimoptions('linprog','TolFun', 10^(-8), 'Display','off')) ; 
[delta_ub,ub] = linprog(-l, A, b, [], [], [], [],  optimoptions('linprog','TolFun', 10^(-8), 'Display','off'));
ub = -ub; %Need to take negative of the minimum to get the max

cs = [lb; ub];

%Compute slack at upper and lower bound solutions
slack_lb = A * delta_lb - b;
slack_ub = A * delta_ub - b;


end
