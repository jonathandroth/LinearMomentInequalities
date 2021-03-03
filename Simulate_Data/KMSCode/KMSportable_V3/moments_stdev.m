function [f_stdev_ineq,f_stdev_eq] = moments_stdev(W,f_ineq,f_eq,J1,J2,KMSoptions)
%% USER-SPECIFIED FUNCTION: Moment Standard Deviation
% The moment functions are in the form
%
%       E_P[m(W_i,theta)] = E_P[f(W_i)] + g(theta)
%
% where
%
%       E_P[m_j(W_i,theta)] <= 0 for j = 1,...,J1
%       E_P[m_j(W_i,theta)] = 0  for j = J1+1,...,J1+J2
%
% This function computes an estimator for the standard deviation of the
% moments: stdev( m(W_i,theta)).  Using separability and noting that theta
% is not random, we have:
%
%        stdev( m(W_i,theta)) = stdev(f(W_i)).
%
% So only f(W_i) is required to compute stdev in the separable case.
%
% INPUT:
%   f_ineq      J1-by-1 vector of moment inequalities f_j(W) for j=1,...,J1
%
%   f_eq        2*J2-by-1 vector of moment inequalities f_j(W) for j=1,...,J2
%               Note that we re-write the moment equalities as two moment
%               inequalities.  Thus, we have f_eq = [f(W) ; - f(W)], where
%               f(W) is the vector of moment equalities.
%
%   J1            Integer number of moment inequalities
%
%   J2.           Integer number of moment equalities
%
%   KMSoptions.   This is a structure of additional inputs.  The user can
%                 add parameters to KMSoptions, say KMSoptions.params,
%                 and call KMSoptions.params in the user-specified
%                 functions.
%                 For example, in the 2x2 entry game, we include the
%                 support for the covariates and the probability that
%                 a particular support point occurs.
%
% OUTPUT:
%   f_stdev_ineq   J1-by-1 vector of stdev(f_j(W)) for the moment
%                  inequalities, with j = 1,...,J1
%
%   f_stdev_eq     2*J2-by-1 vector of stdev(f_j(W)) for the moment
%                  equalities, with j=1,...,J2
%
% Below is a list of examples of moment functions.
% Select DGP
DGP = KMSoptions.DGP;

if DGP == -1
   f_stdev_ineq = std(W).';
   f_stdev_eq   = [];
end


if DGP == 0
   f_stdev_ineq = std(W).';
   f_stdev_eq   = [];
end
if DGP == 1 || DGP == 2 || DGP == 3 || DGP == 4 
    %% Example: Rotated cube (DGP-1 - DGP-4)
    f_stdev_ineq = std(W).';
    f_stdev_eq   = [];

elseif DGP == 5 || DGP == 6 || DGP == 7 || DGP == 8
    %% Example 1: 2-by-2 Entry Game
    f_stdev_ineq(:,1) = sqrt(abs(f_ineq).*(1-abs(f_ineq)));
    f_stdev_eq(:,1)   = sqrt(abs(f_eq).*(1-abs(f_eq)));
end
end
