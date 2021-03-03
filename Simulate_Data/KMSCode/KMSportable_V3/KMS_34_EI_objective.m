function [gamma,dgamma] = KMS_34_EI_objective(theta_aug,KMSoptions)
%% Code description: Expected Improvement with fmincon
% This function computes the objective of the fminimax program using
% fmincon.  The objective function is simply gamma, a constant.  
%
% The objective function is max_{theta,gamma) gamma, and the constraints
% are F_j(theta) <= gamma, where F_j(theta) is the expected improvement of
% the ith moment.
%
% INPUT:
%
% theta_aug     (1+dim_p)-by-1 vector, which includes as its dim_p+1
%               element the parameter gamma.
%
%   KMSoptions.         This is a structure of additional inputs held
%                       constant over the program.  In the 2x2 entry game,
%                       KMSoptions includes the support for the covariates 
%                       and the probability of support point occuring.  
%                       There are also options in KMSoptions to  specify 
%                       optimization algorithm, tolerance, and tuning 
%                       parameters.  However, it is not recommended that 
%                       the user adjusts these.
%  
% OUTPUT:
%   gamma       1-by-1 constant

%   dgamma      (dim_p+1)-by-1 gradient


%% Extract relevant information from KMSoptions
dim_p       = KMSoptions.dim_p;

%% Compute objective and gradient
gamma = theta_aug(end,1);

if nargout > 1
    dgamma = zeros(dim_p+1,1);
    dgamma(end,1) = 1;
end


end