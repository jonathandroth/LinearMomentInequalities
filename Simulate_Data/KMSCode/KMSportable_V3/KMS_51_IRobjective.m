function [val,dval] = KMS_51_IRobjective(theta,q)
%% Code Description
% This function computes the value and gradient of the min/max p'theta
% subject to theta satisfying the population moments.
% INPUTS:
%   theta               dim_p-by-1 parameter vector
%
%   q                   dim_p-by-1 directional vector.  This is either
%                       p or -p
%
% OUTPUT:
%   val                 value of f(theta) = p'theta
%
%   dval                gradient of f(theta) = p'theta


val = -q.'*theta;
ind = find(q~=0);
dval = zeros(size(theta));
dval(ind) = -q(ind);

end