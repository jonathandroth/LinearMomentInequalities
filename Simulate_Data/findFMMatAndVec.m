
function [A_fm,d_fm] = findFMMatAndVec(X)
%This function takes a matrix X and creates matrices A_fm, d_fm
% so that the system: there exists delta s.t. E[Y - Xdelta] \leq 0 is
% equivalent to A_fm * E[Y] \leq d_fm

A = [-X, eye(size(X,1))];
d = zeros(size(X,1),1);
nvar = size(X,1); %number of variables to have at the end, i.e. dim of Y
[A_fm,d_fm] = fourmotz(A,d,nvar);

end