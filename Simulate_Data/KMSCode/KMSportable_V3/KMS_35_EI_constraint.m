function [EI,Ceq,dEI,DCeq] = KMS_35_EI_constraint(theta_aug,q,theta_hash,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,dmodel,KMSoptions)
%% Code description: Expected Improvement with fmincon
% This function computes the objective of the fminimax program using
% fmincon.  The objective function is simply gamma, a constant.  
%
% The objective function is max_{theta,gamma) gamma, and the constraints
% are F_j(theta) <= gamma, where F_j(theta) is the expected improvement of
% the ith moment.
%
% INPUT:
% theta_aug     (dim_p+1)-by-1 parameter vector, which includes gamma in
%               its last component
%
% q             dim_p-by-1 directional vector.  This is either  p or -p.
%
% theta_A       S-by-dim_p matrix of parameter vectors previously explored.    
%
% dmodel        DACE kriging model (computed using theta_A)
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
%   Ei          J-by-1 vector of expected improvements for each moment 
%               inequality minus gamma.  
%
%   dEi         J-by-dim_p matrix of gradients of expected improvement
%               inequality constraints minus 1.

%% Extract relevant information from KMSoptions
J1          = KMSoptions.J1;
J2          = KMSoptions.J2;
J           = KMSoptions.J;     
n           = KMSoptions.n;
dim_p       = KMSoptions.dim_p;

%% Extract theta and gamma
gamma = theta_aug(end,1);
theta = theta_aug(1:dim_p,1);

%% No equality constraints:
Ceq = [];
DCeq = [];

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
h_theta_keep = h_theta;
f_keep = [f_ineq_keep;f_eq_keep];
h_theta_keep(f_keep == 0,:) = [];

% Step 2) c(theta) and s(theta)
% Approximated value of c(theta) using DACE
% If gradient is required, we calcualte gradient of c.
if nargout <= 2
    c_theta    = predictor(theta,dmodel);
else
    [c_theta,dc_theta] = predictor(theta,dmodel);
end

% Compute s^2(theta) 
if nargout <= 2
    [~,~,mse,~]=predictor(theta,dmodel);
    s = sqrt(mse);
else
    [~,~,mse,dmse]=predictor(theta,dmodel);
    s = sqrt(mse);
    ds = 0.5*dmse./s;
end

% Step 3) Compute expected improvement minus gamma
EI = q.'*(theta - theta_hash)*(-normcdf(-(h_theta_keep-c_theta)/s)) - gamma;


%% Gradient (if required)
% The gradient of expected improvement is
%
% dEI_j  = term1 + term2*term3 
%        = (q.').*(-Phi(.)) + q.'(theta - theta#)(-phi(.))* dArg/dtheta)
%
% where dArg/dtheta = -( (dh + dc)s - (h+s)ds)/s^2
%
if nargout > 2
    % First term: (q.').*Phi(.)
    term1 = repmat(q.',[J,1]).*repmat(-normcdf(-(h_theta-c_theta)/s),[1,dim_p]);
    
    % Second term: q.'(theta - theta#)phi(.) 
    term2 = q.'*(theta - theta_hash)*repmat(-normpdf(-(h_theta-c_theta)/s),[1,dim_p]);
    
    % Third term: dArg/dtheta is the derivative of the term inside the
    % argument of the normal CDF.  For this, we need the gradient of the
    % standardized moment conditions, gradient of c(theta), and gradient 
    % MSE.
    
    % Compute gradient of g
    [Dg_ineq,Dg_eq] = moments_gradient(theta,J1,J2,KMSoptions);
    
    % Gradient of Standardized moments:
    Dh_ineq = sqrt(n)*Dg_ineq./repmat(f_stdev_ineq,[1,dim_p]);
    Dh_eq = sqrt(n)*Dg_eq./repmat(f_stdev_eq,[1,dim_p]);
    Dh_theta = [Dh_ineq; Dh_eq];
   
    % Derivative of the term inside the arguement of the CDF:
    term3 = -((Dh_theta - repmat(dc_theta.',[J,1]))*s - repmat(h_theta-c_theta,[1,dim_p]).*repmat(ds.',[J,1]))/s^2;
    
    % Gradient of EI
    DEI = term1 + term2.*term3;
    
    % Drop moments close to boundary
    DEI(f_keep == 0,:) = [];
    
    % Transpose to get correct orientation for fmincon
    DEI = DEI.';
    
    % Include row of -1 to get correct gradient of constraint function
    num_keep = size(find(f_keep == 1),1);
    dEI = [DEI; -1*ones(1,num_keep)];
end


end