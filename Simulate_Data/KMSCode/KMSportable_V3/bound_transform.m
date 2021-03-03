function [LB_out,UB_out] = bound_transform(LB_in,UB_in,KMSoptions)
%% USER-SPECIFIED FUNCTION
% This function is only relevant when polytope constraints are included and
% the draw-and-discard sampling method is used (HR = 0).  In the
% contraction step, and can become difficult to find points that satisfy
% the polytope constraints.  In certain circumstances we use information
% about the polytope constraints to reduce the size of the box constraints
% and make it more likely that our points satisfy the polytope constraints.
%
% Example: DGP8.  The polytope constraints in DGP8 require that 
%   theta_j <= min(theta_1,theta_2) for j=3,4,5.
% So we can transform UB_out by setting UB_out to the minimum of UB_in and
% min(theta_1,theta_2).
%
% INPUT:
%   LB_in               Lower bound input
%   UB_in               Upper bound input
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
%   LB_out              Transformed lower bound 
%   UB_out              Transformed upper bound 

% DGP8 example
DGP = KMSoptions.DGP;
UB_out = zeros(size(UB_in));
LB_out = zeros(size(LB_in));

% LB_out is the same as LB_in
LB_out = LB_in;

if DGP == 8
    min_UB       = min(UB_in(1:2));
    UB_out(3:5)  = min( UB_in(3:5), min_UB);
    UB_out(1:2)  = UB_in(1:2);
    UB_out       = max(LB_out,UB_out);          % Make sure UB >= LB.
else
    UB_out = UB_in;
end



end