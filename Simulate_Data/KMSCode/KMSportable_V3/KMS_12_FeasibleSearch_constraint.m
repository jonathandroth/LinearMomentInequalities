function [m_theta,Ceq,Dm_theta,DCeq] = KMS_12_FeasibleSearch_constraint(theta_aug,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,KMSoptions)
%% Code description: Expected Improvement with fmincon
% This function computes the objective of the fminimax program using
% fmincon.  The objective function is simply gamma, a constant.  
%
% The objective function is max_{theta,gamma) gamma, and the constraints
% are F_j(theta) <= gamma, where F_j(theta) is the expected improvement of
% the ith moment.
%
% INPUT:
%   theta               dim_p-by-1 parameter vector
%                       
%   f_ineq,f_eq         Empirical moments
%
%   f_ineq_keep,f_eq_keep       Moments to keep
%
%   f_stdev_ineq,f_stdev_eq     Standard deviation of empirical moments
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
%   m_theta             J-by-1 vector of standardized moment inequalities
%                       minus gamma (control parameter)
%
%   Dm_theta            dim_p-by-J matrix of gradients of standardized
%                       moment inequalities   

%% Extract relevant information from KMSoptions
J1          = KMSoptions.J1;
J2          = KMSoptions.J2;
J           = KMSoptions.J;     
n           = KMSoptions.n;
dim_p       = KMSoptions.dim_p;

%% Extract theta and gamma
gamma = theta_aug(end,1);
theta = theta_aug(1:dim_p,1);

%% No inequality constraints:
Ceq = [];
DCeq = [];

%% Standardized moments
% Theoretical moments:
[g_ineq,g_eq] = moments_theta(theta,J1,J2,KMSoptions);

% Standardized momoments using empirical and theoretical moments:
m_theta = sqrt(n)*(([f_ineq;f_eq] + [g_ineq;g_eq])./[f_stdev_ineq;f_stdev_eq]);

% Drop moments close to boundary
f_keep = [f_ineq_keep;f_eq_keep];
m_theta(f_keep == 0,:) = [];

% Step 3) Compute standardized moments minus gamma
m_theta = m_theta - gamma;


%% Gradient (if required)
% Output gradient iff it is requested

if nargout > 2
    % Theoretical moment gradients
    [Dg_ineq,Dg_eq] = moments_gradient(theta,J1,J2,KMSoptions);
    Dg = [Dg_ineq ; Dg_eq];
    
    % Standardized moment gradient:
    Dm_theta = sqrt(n).*Dg./repmat([f_stdev_ineq;f_stdev_eq],1,dim_p);
    
    % Drop moments close to binding
    Dm_theta(f_keep == 0,:) = [];
     
    % Transpose for fminimax  
    Dm_theta = Dm_theta.';
    
    % Include row of -1 to get correct gradient of constraint function
    num_keep = size(find(f_keep == 1),1);
    Dm_theta = [Dm_theta; -1*ones(1,num_keep)];
end


end