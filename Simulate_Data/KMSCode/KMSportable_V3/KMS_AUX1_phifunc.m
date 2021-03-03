function [ value ] = KMS_AUX1_phifunc(x)
%% Code Description: GMS hard thresholding function, phi
% This function computes the hard thresholding GMS function
%
%       phi(x) = 0 if x >= -1 and phi(x) = -infinity if x < -1.
%
% INPUT:
%       x       J-by-1 vector, where J is total number of moments
%
% OUTPUT:
%   value       J-by-1 vector, where value(jj) = 0 iff x(jj) >= -1 and
%               value(jj) = -infinity if x(jj) < -1.

value = zeros(size(x));
value(x<-1) = -inf;

end

