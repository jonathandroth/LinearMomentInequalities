function [EI] = KMS_37_EI_value(theta,q,theta_hash,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,dmodel,KMSoptions)
%% Code description: Expected Improvement with fmincon
% This function computes the expected improvement at theta.
%
% INPUT:
% theta         (dim_p)-by-1 parameter vector
%
% q             dim_p-by-1 directional vector.  This is either  p or -p.
%
% dmodel        DACE kriging model (computed using theta_A)
%
% KMSoptions    This is a structure of additional inputs held 
%               constant over the program.  In the 2x2 entry game, 
%               KMSoptions includes the support for the covariates
%               and the probability of support point occuring. 
%  
% OUTPUT:
%   Ei          J-by-1 vector of expected improvements for each moment 
%               inequality minus gamma.  

%% Extract relevant information from KMSoptions
J1          = KMSoptions.J1;
J2          = KMSoptions.J2;
n           = KMSoptions.n;


%% Expected Improvement
% We compute expected improvement for each moment inequality j=1...J:
%   EI_j = (q'theta - q'theta_#)_{+}*Phi( (h_j(theta) - c(theta))/s(theta))
% where h_j(theta) is the standardized moment
%   h_j(theta) = sqrt(n)*m_j(X,theta)/sigma(X)
% and c_L(theta), s_L(theta) are from the DACE auxillary model.
% Note that we are searching over the space of theta such that 
%   q'theta >= q'theta_#
% so the max(0,.) is not required.

% Step 1) h_j(theta)
% We compute the standardized moments
% Theoretical momoments
[g_ineq,g_eq] = moments_theta(theta,J1,J2,KMSoptions);

% Standardized momoments
h_theta = sqrt(n)*(([f_ineq;f_eq] + [g_ineq;g_eq])./[f_stdev_ineq;f_stdev_eq]);

% Drop moments close to boundary
f_keep = [f_ineq_keep;f_eq_keep];
h_theta(f_keep == 0,:) = [];

% Step 2) c(theta) and s(theta)
% Approximated value of c(theta) using DACE
% If gradient is required, we calcualte gradient of c.
c_theta    = predictor(theta,dmodel);

% Compute s^2(theta) 
[~,~,mse,~]=predictor(theta,dmodel);
s = sqrt(mse);

% Step 3) Compute expected improvement minus gamma
EI = q.'*(theta - theta_hash)*(-normcdf(-(h_theta-c_theta)/s));

end