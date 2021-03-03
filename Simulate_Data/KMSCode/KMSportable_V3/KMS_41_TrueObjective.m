function [val,dval] = KMS_41_TrueObjective(theta,q)
%% Code Description
% This function computes the value and gradient of the min/max problem on
% Pg 13, Eq 2.16.   We use fmincon to find the upper/lower bounds of the
% confidence set.  This method is likely to perform poorly in realistic
% problems.  This method is used as an error-check in simple problems 
% to ensure that EAM is working correctly.

% INPUTS:
%   theta               dim_p-by-1 parameter vector
%
%   q                   dim_p-by-1 directional vector.  This is either
%                       p or -p
%
% OUTPUT:
%   val                 value of problem 2.16: p'theta
%
%   dval                gradient of problem 2.16


val = -q.'*theta;
ind = find(q~=0);
dval = zeros(size(theta));
dval(ind) = -q(ind);

end