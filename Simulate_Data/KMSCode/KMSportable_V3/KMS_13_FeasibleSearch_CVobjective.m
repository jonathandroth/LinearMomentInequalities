function [val,Dval] = KMS_13_FeasibleSearch_CVobjective(theta,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,KMSoptions)
%% Code description: Expected Improvement with fmincon
% The objective function outputs the constraint violation 
% determined by the sum of differences in the constraints.  Minimizing the
% constraint violation has proven to be unreliable, so it is currently not
% in use.  In principle, we can try find a feasible point using both this
% method and the fminimax method.
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
%   val             Value of objective function, which is the sum over j of
%                   max(0, normalized moment)^2.
%
%   Dval            dim_p-by-J matrix of gradients of standardized
%                       moment inequalities   

%% Extract relevant information from KMSoptions
J1          = KMSoptions.J1;
J2          = KMSoptions.J2;
J           = KMSoptions.J;     
n           = KMSoptions.n;
dim_p       = KMSoptions.dim_p;

%% Standardized moments
% Theoretical moments:
[g_ineq,g_eq] = moments_theta(theta,J1,J2,KMSoptions);

% Standardized momoments using empirical and theoretical moments:
m_theta = sqrt(n)*(([f_ineq;f_eq] + [g_ineq;g_eq])./[f_stdev_ineq;f_stdev_eq]);

% Drop moments close to boundary
f_keep = [f_ineq_keep;f_eq_keep];
m_theta(f_keep == 0,:) = [];

% Step 3) Compute value of the objective function
val = sum(max(0,m_theta).^2);


%% Gradient (if required)
% Output gradient iff it is requested

if nargout > 1
    % Theoretical moment gradients
    [Dg_ineq,Dg_eq] = moments_gradient(theta,J1,J2,KMSoptions);
    Dg = [Dg_ineq ; Dg_eq];
    
    % Standardized moment gradient:
    Dm_theta = sqrt(n).*Dg./repmat([f_stdev_ineq;f_stdev_eq],1,dim_p);
    
    % Drop moments close to binding
    Dm_theta(f_keep == 0,:) = [];
    
    % Transpose for fminimax 
    %Dm_theta = Dm_theta.';
    % Dm_theta is now dim_p-by-J.
    
    % Gradient is equal to sum_j 2*max(0, m_j(theta))*Dm_j(theta).
    Dval = sum(repmat(2*max(0,m_theta),[1,dim_p]).*Dm_theta).';
    
end


end