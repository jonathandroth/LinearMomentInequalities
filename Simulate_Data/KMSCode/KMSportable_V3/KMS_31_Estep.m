function [c_Estep,CV_Estep,max_violation] = KMS_31_Estep(theta_Estep,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions)
%% Code Description: Evaluation Step
% This function executes the E-step in the EAM algorithm.  See Pg 10-12,
% and in particular:
%       Pg 10, eq 2.10-2.11
%       Pg 11, fixed-point algorithm (Brent's Method is used)
%       Pg 12, description of E-step
%
% This function takes in as an input, among other things, theta_Estep.
% theta_Estep is a dim_e-by-dim_p matrix, where each row represents a
% vector theta in the parameter space [theta_lo,theta_hi].  For each
% theta_l, we compute the (1-alpha)-critical value (see eq 2.10 and
% 2.11). Operationaly, we compute the critical value using the algorithm 
% Pg 11.  
%
% NOTE: The code now uses Brent's fixed-point algorithm for accelerated
% convergence.
%
% INPUTS:
%   theta_Estep         K-by-dim_p matrix.  Each row is a parameter vector
%                       from the parameter space [theta_lo,theta_hi].
%
%   f_ineq, f_eq        Empirical moments
%
%   f_ineq_keep,f_eq_keep       Moments to keep
%
%   f_stdev_ineq,f_stdev_eq     Standard deviation of moments
%
%   G_ineq,G_eq         Bootstrapped and recentered moments
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

%   c_Estep             dim_e-by-1 vector of(1-alpha) critical value 
%                       evaluated at theta_l, l = 1,...,dim_e.    
%
%   CV_Estep            dim_e-by-1 vector that determines the constraint
%                       violation.  theta_l is feasible iff CV_Estep =0.
%
%   max_violation       Maximum violation of c -standardized moments

%% Extract relevant information from KMSoptions
parallel    = KMSoptions.parallel;  
J1          = KMSoptions.J1;
J2          = KMSoptions.J2;
J3          = KMSoptions.J3;
J           = KMSoptions.J;
kappa       = KMSoptions.kappa;
n           = KMSoptions.n;
phi         = KMSoptions.phi; 
paired_mom  = KMSoptions.paired_mom;
dim_p       = KMSoptions.dim_p; 

%% Preliminary
dim_e =  size(theta_Estep,1);        % Number of evaluation points
c_Estep = zeros(dim_e,1);            % Critical value for each evaluation point
CV_Estep = zeros(dim_e,1);
max_violation = zeros(dim_e,1);

%% Calculate critical value for each theta_l, l=1,...,dim_e
if parallel == 1
    % Run parfor if we are running parallel program
    parfor ll = 1:dim_e
        % 0) Select theta_l
        % Fix a test value of theta, theta_l.
        % theta_l is passed through as a column vector.
        theta_test = theta_Estep(ll,:).';
        
        % 1) GMS function
        % (See Pg 9, Section 2.2, Eq 2.8-2.9)
        % The GMS function, phi, is a J-by-1 vector that is (for the case
        % of hard thesholing) equal to 0 iff moment j is close to binding
        % and equal to -infinity iff moment j is far from binding.  All
        % moment equalities are binding, so phi_j = 0 for the moment
        % equalities. To determine whether or not a moment is close to
        % binding, we use the function xi_j given in Eq 2.8.
        phi_test = zeros(J,1);
        
        % Compute theoretical bounds g(theta).
        % These theoretical enter the GMS function in Eq 2.8.  
        [g_ineq,g_eq] = moments_theta(theta_test,J1,J2,KMSoptions);
        
        % Measure of close to binding (Equation 2.8)
        % Note: moment is m(W,theta) = f(W) + g(theta)
        % So measure of "close to binding" is given by
        % xi = (1/kappa)*sqrt(n)*(f(W) + g(theta))/(stdev)
        % for the moment inequalities, and zero for the moment equalities
        xi_ineq = (1/kappa)*sqrt(n).*(f_ineq + g_ineq)./f_stdev_ineq;
        xi = [xi_ineq; zeros(2*J2,1)];

        % GMS function (Equation 2.9)
        % phi_test is the GMS function evaluated at the measure of "close
        % to binding" xi, computed above.  
        % In the KMS paper, the hard-threshing GMS function is used,
        % where phi = -inifity if xi < -1, 0 else.
        phi_test = phi(xi);
        
        % 2) Gradients
        % (See Pg 8, Eq 2.6)
        % We compute the gradient of (m(W,theta)/sig(W,theta)).  
        % Recall m(W,theta)/std(W) = (f(W) + g(theta))/std(W).  
        % So the gradient of m(W,theta)/std(W) is equal to the gradient of
        % g(theta), Dg(theta) divided by the standard error. 
        [Dg_ineq ,  Dg_eq] = moments_gradient(theta_test,J1,J2,KMSoptions);

        % Gradients are normalized by the standard error
        Dg_ineq = Dg_ineq./repmat(f_stdev_ineq,[1,dim_p]);
        Dg_eq = Dg_eq./repmat(f_stdev_eq,[1,dim_p]);

        % 3) Paired moment inequalities
        % The paired moments are those that are highly correlated.  We
        % assume that the only moment inequalities that are highly
        % correlated come in pairs.  That is, we do not allow for moments
        % 1,2, and 3 to be highly correlated. 
        % If two paired moment inequalities are "close to binding", then 
        % we force both the moment and the gradient to be equal to the
        % negative of one another (i.e., we force the moment inequalities
        % to become one moment equality).  Note that, however, in the 
        % separable case this only impacts the gradients rather than the
        % moments, since the Gaussian recentered process does not depend on
        % theta.  
        for jj = 1:J3
            % Find the jjth paired moment inequality:
            ind = find(paired_mom == jj); 

            % Check to see if both are selected by GMS
            if phi_test(ind(1)) == 0 && phi_test(ind(2)) == 0
               % In this case, both are selected by GMS  so we adjust
               % moments and gradients.
               % In the separable case we only need to adjust gradients and not
               % moments.
               Dg_ineq(ind(2),:) = -Dg_ineq(ind(1),:);
            end
        end
        
        % 4) Compute rho polytope constraints
        [A_rho,b_rho] = Rho_Polytope_Box(theta_test,KMSoptions);
        
        % 5) Compute critical value
        % This follows the algorithm on page 11. Essentially, the 
        % algorithm is a root-finding algorithm that finds the root to the
        % function h(c) = (1/n) sum_b psi_b(c) - (1-alpha).  If h(c) = 0,
        % then the coverage of 1-alpha is obtained at theta_test.
        c_Estep(ll,1) = KMS_32_Critval(phi_test,f_ineq_keep,f_eq_keep,G_ineq,G_eq,Dg_ineq,Dg_eq,A_rho,b_rho,KMSoptions);
        
        % 6) Constraint violation
        % Standardized moments
        m_theta = sqrt(n)*(([f_ineq;f_eq] + [g_ineq;g_eq])./[f_stdev_ineq;f_stdev_eq]);
        
        % Drop moments with value of f(W) close to boundary
        f_keep = [f_ineq_keep;f_eq_keep];
        m_theta(f_keep == 0,:) = [];
    
        % Constraint violation:
        CV_Estep(ll,1) = sum(max(0,m_theta-c_Estep(ll,1)).^2);
        
        % Maximum violation:
        max_violation(ll,1) = max(m_theta - c_Estep(ll,1));
    end
else
    % If parallel computing not requested, run for loop.
    for ll = 1:dim_e
        % 0) Select theta_l
        % Fix a test value of theta, theta_l.
        % theta_l is passed through as a column vector.
        theta_test = theta_Estep(ll,:).';
        
        % 1) GMS function
        % (See Pg 9, Section 2.2, Eq 2.8-2.9)
        % The GMS function, phi, is a J-by-1 vector that is (for the case
        % of hard thesholing) equal to 0 iff moment j is close to binding
        % and equal to -infinity iff moment j is far from binding.  All
        % moment equalities are binding, so phi_j = 0 for the moment
        % equalities. To determine whether or not a moment is close to
        % binding, we use the function xi_j given in Eq 2.8.
        phi_test = zeros(J,1);
        
        % Compute theoretical bounds g(theta).
        % These theoretical enter the GMS function in Eq 2.8.  
        [g_ineq,g_eq] = moments_theta(theta_test,J1,J2,KMSoptions);
        
        % Measure of close to binding (Equation 2.8)
        % Note: moment is m(W,theta) = f(W) + g(theta)
        % So measure of "close to binding" is given by
        % xi = (1/kappa)*sqrt(n)*(f(W) + g(theta))/(stdev)
        % for the moment inequalities, and zero for the moment equalities
        xi_ineq = (1/kappa)*sqrt(n).*(f_ineq + g_ineq)./f_stdev_ineq;
        xi = [xi_ineq; zeros(2*J2,1)];

        % GMS function (Equation 2.9)
        % phi_test is the GMS function evaluated at the measure of "close
        % to binding" xi, computed above.  
        % In the KMS paper, the hard-threshing GMS function is used,
        % where phi = -inifity if xi < -1, 0 else.
        phi_test = phi(xi);
        
        % 2) Gradients
        % (See Pg 8, Eq 2.6)
        % We compute the gradient of (m(W,theta)/sig(W,theta)).  
        % Recall m(W,theta)/std(W) = (f(W) + g(theta))/std(W).  
        % So the gradient of m(W,theta)/std(W) is equal to the gradient of
        % g(theta), Dg(theta) divided by the standard error. 
        [Dg_ineq ,  Dg_eq] = moments_gradient(theta_test,J1,J2,KMSoptions);

        % Gradients are normalized by the standard error
        Dg_ineq = Dg_ineq./repmat(f_stdev_ineq,[1,dim_p]);
        Dg_eq = Dg_eq./repmat(f_stdev_eq,[1,dim_p]);

        % 3) Paired moment inequalities
        % The paired moments are those that are highly correlated.  We
        % assume that the only moment inequalities that are highly
        % correlated come in pairs.  That is, we do not allow for moments
        % 1,2, and 3 to be highly correlated. 
        % If two paired moment inequalities are "close to binding", then 
        % we force both the moment and the gradient to be equal to the
        % negative of one another (i.e., we force the moment inequalities
        % to become one moment equality).  Note that, however, in the 
        % separable case this only impacts the gradients rather than the
        % moments, since the Gaussian recentered process does not depend on
        % theta.  
        for jj = 1:J3
            % Find the jjth paired moment inequality:
            ind = find(paired_mom == jj); 

            % Check to see if both are selected by GMS
            if phi_test(ind(1)) == 0 && phi_test(ind(2)) == 0
               % In this case, both are selected by GMS  so we adjust
               % moments and gradients.
               % In the separable case we only need to adjust gradients and not
               % moments.
               Dg_ineq(ind(2),:) = -Dg_ineq(ind(1),:);
            end
        end

        % 4) Compute rho polytope constraints
        [A_rho,b_rho] = Rho_Polytope_Box(theta_test,KMSoptions);

        % 5) Compute critical value
        % This follows the algorithm on page 11. Essentially, the 
        % algorithm is a root-finding algorithm that finds the root to the
        % function h(c) = (1/n) sum_b psi_b(c) - (1-alpha).  If h(c) = 0,
        % then the coverage of 1-alpha is obtained at theta_test.
        c_Estep(ll,1) = KMS_32_Critval(phi_test,f_ineq_keep,f_eq_keep,G_ineq,G_eq,Dg_ineq,Dg_eq,A_rho,b_rho,KMSoptions);

        % 6) Constraint violation
        % Standardized moments
        m_theta = sqrt(n)*(([f_ineq;f_eq] + [g_ineq;g_eq])./[f_stdev_ineq;f_stdev_eq]);
        
        % Drop moments close to boundary
        f_keep = [f_ineq_keep;f_eq_keep];
        m_theta(f_keep == 0,:) = [];
    
        % Constraint violation:
        CV_Estep(ll,1) = sum(max(0,m_theta-c_Estep(ll,1)).^2);
        
        % Maximum violation:
        max_violation(ll,1) = max(m_theta - c_Estep(ll,1));
    end
end
end

