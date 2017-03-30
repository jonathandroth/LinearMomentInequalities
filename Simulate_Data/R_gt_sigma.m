
%This function takes as input
% g_T: a vector of moments
% Sigma: a covariance matrix

%It returns:

% R( g_T, Sigma), as defined in Isaiah's "Confidence Sets..." Notes

function R = R_gt_sigma( g_T, Sigma)

R = max( g_T ./ sqrt( diag( Sigma) ) ); 

end
