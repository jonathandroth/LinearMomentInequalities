function [A,b] = Rho_Polytope_Box(theta,KMSoptions)
%% USER-SPECIFIED FUNCTION: Rho Polytope Constraints
% This is a user-written function that may be necessary if the user
% specifies polytope constraints rather than simple box constraints.  We
% require this function for DGP8, since at certain values of theta, the
% rho-box around theta may intersect the polytope constraint.  This
% function prevents this issue by taking a rho-polytope around theta. 
%
% S = number of additional polytope constraints
%
% Please set A = [], b = [] if there are is no user-specified rho-polytope
% constraints.
%
% INPUT:
%   theta       dim_p-by-1 parameter vector 
%
%   KMSoptions  This is a structure of additional inputs.  The user can
%               add parameters to KMSoptions, say KMSoptions.params,
%               and call KMSoptions.params in the user-specified
%               functions.
%               For example, in the 2x2 entry game, we include the
%               support for the covariates and the probability that
%               a particular support point occurs.
% OUTPUT:
%   A           S-by-dim_p polytope constraint
%
%   b           S-by-1 polytope constraint

DGP     =   KMSoptions.DGP;


if DGP == 8
    % See DGP8 note.
    % We include constraints that prevent the rho-box from hitting the
    % boundary on the polytope parameter space
    
    dim_p   =   KMSoptions.dim_p;
    S       =   KMSoptions.S;
    n       =   KMSoptions.n;
    A = zeros(S,dim_p);
    b = zeros(S,1);
    
    A(1:2 ,1:2)  = eye(2);
    A(3:4 ,1:2)  = -eye(2);
    A(5:7 ,3:5)  = -eye(3);
    A(8:10,1)    = -ones(3,1);
    A(8:10,3:5)  = eye(3);
    A(11:13,2)   = -ones(3,1);
    A(11:13,3:5) = eye(3);
    
    b = [sqrt(n)*(1-theta(1))   ;
        sqrt(n)*(1-theta(2))    ;
        sqrt(n)*theta(1)        ;
        sqrt(n)*theta(2)        ;
        sqrt(n)*theta(3)        ;
        sqrt(n)*theta(4)        ;
        sqrt(n)*theta(5)        ;
        -sqrt(n)*(-theta(1) + theta(3))  ;
        -sqrt(n)*(-theta(1) + theta(4))  ;
        -sqrt(n)*(-theta(1) + theta(5))  ;
        -sqrt(n)*(-theta(2) + theta(3))  ;
        -sqrt(n)*(-theta(2) + theta(4))  ;
        -sqrt(n)*(-theta(2) + theta(5))  ];
    
else
   A = [];
   b = [];
end

