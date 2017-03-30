% This function takes a matrix of size T_m x M, where T_m is the number of
% independent moments and M is the number of moments. Each column is
% a vector with the independent observations of the given moment

%The function returns g_t and Sigma as in Isaiah's note "Confidence Sets
%for Parameters..." 
% g_T:
% Sigma

function [g_T, Sigma] = gT_and_Sigma_fn( moment_matrix )

T_m = size(moment_matrix,1);

g_T = 1/sqrt(T_m) *  sum( moment_matrix, 1)' ;
Sigma = cov(moment_matrix); 

end
