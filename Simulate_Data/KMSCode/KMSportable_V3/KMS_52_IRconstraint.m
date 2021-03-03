function [C,Ceq,dC,dCeq] = KMS_52_IRconstraint(theta,f_ineq_pop,f_eq_pop,KMSoptions)
% This function computes the constraint of the min/max problem p'theta
% subject to theta satisfying the population moments.
%
% INPUTS:
%   theta               dim_p-by-1 parameter vector
%
%   f_ineq_pop,f_eq_pop  population moments
%
%   KMSoptions.         This is a structure of additional inputs held
%                       constant over the program.  In the 2x2 entry game,
%                       KMSoptions includes the support for the covariates
%                       and the probability of support point occuring.
%
% OUTPUT:
%   c                   J-by-1 vector of inequality constraints
%
%   dc                  Gradient

%% No equality constraints


%% Extract relevant information from KMSoptions
J1          = KMSoptions.J1;
J2          = KMSoptions.J2;  

%% Moments
[g_ineq,g_eq] = moments_theta(theta,J1,J2,KMSoptions);

% Moments
C= [f_ineq_pop] + [g_ineq];
Ceq = f_eq_pop(1:J2,:) + g_eq(1:J2,:);


%% Gradient if requested
if nargout>2
    [Dg_ineq,Dg_eq] = moments_gradient(theta,J1,J2,KMSoptions);
    dC = [Dg_ineq];
    dC = dC.';
    dCeq = Dg_eq(1:J2,:);
    dCeq = dCeq.';
end


end