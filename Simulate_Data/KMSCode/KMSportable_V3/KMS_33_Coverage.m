function [psi_1] = KMS_33_Coverage(c_1,psi_0,psi_2,phi_test,f_ineq_keep,f_eq_keep,G_ineq,G_eq,Dg_ineq,Dg_eq,A_rho,b_rho,KMSoptions)
%% Coverage Function
% This function computes the psi function in Equation (2.12) evaluated at 
% c_1.  psi(c) is a Bx1 vector of indicator functions.  The bbth element in
% psi(c) is equal to one iff the constraint set intersects the hyperplane
% p'lambda = 0.
%
% INPUTS:
%   c_1                 Current iteration's critical value
%
%   psi_0               Previous iteration's vector of indicators: psi(c_0)
%
%   psi_2               Previous iteration's vector of indicators: psi(c_2)
%
%   phi_test            GMS function evaluated at theta_test
% 
%   f_ineq_keep,f_eq_keep   Moments to keep
%
%   G_ineq,G_eq         Bootstrapped and recentered moments
%
%   Dg_ineq,Dg_eq       Gradients of moment function evaluated at theta_test
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
%   psi_1               Bx1 vector of indicators.  The bbth element is an
%                       indicator equal to one iff the constraint set
%                       evaluated at c_1 intersects the hyperplane 
%                       p'lambda=0
%                       
% Note: The vector psi_1 is updated based c_1 and c_0 for computational
% reasons.  One could compute psi_1 using just c_1, but that is slower.
%    Case 1:  If c_1 > c_0, then psi_1 >= psi_0.
%             So if psi_0=1, then we do not need to recompute
%             psi_1 -- it is equal to 1.  Therefore we only
%             update psi_1 such that psi_0 = 0
%   Case 2:   If c_1 < c_0, then psi_1 <= psi_0.
%             So if psi_0=0, then we do not need to recompute
%             psi_1 -- it is equal to 0.  Therefore we only
%             update psi_1 such that psi_0 = 1

%% Extract relevant information from KMSoptions
rho         = KMSoptions.rho;
component   = KMSoptions.component;
CVXGEN      = KMSoptions.CVXGEN;
type        = KMSoptions.type;
dim_p       = KMSoptions.dim_p;
J           = KMSoptions.J;
CVX_maxit   = KMSoptions.CVX_maxit;
CVXGEN_name = KMSoptions.CVXGEN_name;
S           = KMSoptions.S;

%% Preliminary
p = zeros(dim_p,1);
p(component) = 1;
psi_1 = -999*ones(size(psi_0));  
G = [G_ineq;G_eq];
D = [Dg_ineq;Dg_eq];

% Function handle to call CVXGEN
if ~isempty(CVXGEN_name)
    call_cvx = strcat('@(params,settings)',CVXGEN_name,'(params,settings)');
    call_cvx = str2func(call_cvx);
end

%% Update psi_1 based on psi_0,psi_2
% The assumption is that c_0 < c_1 < c_2.  By monotonicty of psi, it
% follows that psi_0 <= psi_1 <= psi_2.  Since psi is a binary-valued
% vector, we yield the following observations:
%   1) If psi_0 = psi_2 = 0, then psi_1 = 0
%   2) If psi_0 = psi_2 = 1, then psi_1 = 1
% Therefore, it remains to find psi_1(bb) such that psi_0(bb) = 0,
% psi_2(bb) = 1.
psi_1( psi_0 == psi_2 & psi_0 == 0,1) = 0;
psi_1( psi_0 == psi_2 & psi_0 == 1,1) = 1;

% Remaining entries to update:
ind_psi = find(psi_1 == -999);
d_ind = size(ind_psi,1);

%% Run CVXGEN or CVX
% We run CVXGEN or CVX to compute whether or not the constraint set
% intersects p'lambda=0.

% Create parameters to pass to CVXGEN
A = zeros(J + 2*dim_p + 2 + S,dim_p);
b = zeros(J + 2*dim_p + 2 + S,1);

% First J entries of (A,b) are the G + Dx + phi_test <= c or Dx <= c -
% G - phi_test
% CVXGEN cannot handle -inf, so if phi_test = -inf, we set that
% constraint to 0 <= 0
% Nb: We allow for phi_test to be not necessarily equal to 0 (e.g., if
% not using hard-thresholding), so it is important to include this in
% the constraint set.
% 
% We also drop moments that have f(W) close to boundary.
f_keep = [f_ineq_keep;f_eq_keep];
ind_phi = find(phi_test > -inf & f_keep == 1);
A(ind_phi,:) = D(ind_phi,:);
%b(ind,:) = c_1-G(ind,ii)-phi_test(ind,1); 
% NB: b is updated in the loop below, since it depends on bootstrap ii

% Next 2*dim_p entries of (A,b) require that lambda in in a rho-box, so
% that lambda <= rho and -rho <= lambda (i.e., -lambda <= rho).
A(J+1:J+dim_p,:) = eye(dim_p,dim_p);
A(J+dim_p + 1:J+2*dim_p,:) = -eye(dim_p,dim_p);
b(J+1:J+dim_p,1) = rho;
b(J+dim_p + 1:J+2*dim_p,1) = rho;

% Add in rho polytope constraints
A(J + 2*dim_p + 2 + 1 : J + 2*dim_p + 2 + S, :) = A_rho;
b(J + 2*dim_p + 2 + 1 : J + 2*dim_p + 2 + S, 1) = b_rho;

% Finally, the last 2 entries of (A,b) require that p'lambda = 0 in the
% case of two-sided testing.  Otherwise, we have -p'lambda <= 0 for
% one-sided-UB test and p'lambda <= 0 for one-sided-LB test.
if strcmp('two-sided',type) == 1
    A(J+2*dim_p+1,:) = -p;
    A(J+2*dim_p+2,:) = p;
elseif strcmp('one-sided-UB',type) == 1
    A(J+2*dim_p+1,:) = -p;
elseif strcmp('one-sided-LB',type) == 1
    A(J+2*dim_p+2,:) = p;
end

for ii = 1:d_ind
    % Update b (Gaussian process depends on bootstrap ii)
    b(ind_phi,1) = c_1-G(ind_phi,ind_psi(ii))-phi_test(ind_phi,1); 

    % Add Gaussian from bootstrap bb to the structure
    params.A = A;
    params.b = b;
    settings.verbose = 0;                                                     
    settings.eps = 1e-9;   
    settings.max_iters = CVX_maxit;
   
    if CVXGEN == 1
        % Pass structure of parameters to CVXGEN
        [~, status] = call_cvx(params,settings);
    else
        % Pass to CVX
        % In order to speed up convergence, we drop rows of A,b that either 
        % correspond to:
        %   - Moment j is close to binding
        %   - Hardthresholding GMS function did not select moment j
        % We cannot do this in CVXGEN, as we need to preset the
        % number of moments.  
        ind_drop = find(phi_test == -inf | f_keep == 0);
        A_cvx = A;
        b_cvx = b;
        A_cvx(ind_drop,:) = [];
        b_cvx(ind_drop,:) = [];
        
        cvx_begin quiet
            variable lambda(dim_p, 1);
            minimize(0);
            subject to
            A_cvx*lambda <= b_cvx;
        cvx_end
        status.converged    = strcmp(cvx_status, 'Solved');
    end
    
    % If the constraint set intersects p'lambda = 0, then 
    % status.converged
    if status.converged==1  
        psi_1(ind_psi(ii)) = 1;
    else
        psi_1(ind_psi(ii)) = 0;
    end
end

end












