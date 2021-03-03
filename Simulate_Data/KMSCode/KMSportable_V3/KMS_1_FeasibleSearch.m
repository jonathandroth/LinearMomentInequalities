function [theta_feas,flag_feas] =  KMS_1_FeasibleSearch(p,theta_0,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions)
%% Code description: EAM Feasible Search
%  This function executes an auxiliary search for a feasible point.
%  The algorithm attempts to find a feasible point by minimizing the 
%  "hard moment" constraint:
%
%    min max_{j=1,...,J} sqrt(n)m(X,theta)/sigma(X)
%
% In the case that the optimal point from this search satisfies 
%
%    CV = sum_j max(0, (sqrt(n)m(X,theta)/sigma(X) - c(theta))^2) = 0
% 
% we take theta_feas to be this feasible point.  If no such point exists,
% we output empty set, and the KMS confidence interval will be the
% parameter space.  
%
% INPUTS:
%   p                   dim_p-by-1 directional vector.
%
%   theta_0             Initial guess for true theta
%
%   f_ineq,f_eq         Empirical moments
%
%   f_ineq_keep,f_eq_keep       Moments to keep
%
%   f_stdev_ineq,f_stdev_eq     Standard deviation of empirical moments
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
%   theta_feas          dim_p-by-2 parameter vector of feasible points.   
%                       The first row is the point that maximizes p'theta 
%                       among the set of feasible points that we find, and
%                       the second row is the minimizer.
%
%   flag_feas           equal to one iff we find a feasible point.

%% Extract relevant information from KMSoptions
LB_theta           = KMSoptions.LB_theta;
UB_theta           = KMSoptions.UB_theta;
A_theta            = KMSoptions.A_theta;
b_theta            = KMSoptions.b_theta;
sample_method      = KMSoptions.sample_method;
dim_p              = KMSoptions.dim_p;
options_fmincon    = KMSoptions.options_fmincon;
multistart_num     = KMSoptions.multistart_num;
parallel           = KMSoptions.parallel;

%% Feasible Search 
% The fminimax search find points that minimize the "hard" moment 
% conditions:
%
%       min max_{j=1,...,J} sqrt(n)m(X,theta)/sigma(X).
%
% The constraint evaluation is calculated for the solutions to to above
% problem:
%
%       CV = sum_j max(0, [sqrt(n)m(X,theta)/sigma(X) - c(theta)]^2)
%
% We use as initial points in EAM those that satisfy CV = 0.  

% Draw points from Theta and in a parfor loop pass these to fmincon.
if sample_method == 0 || sample_method == 1 
    % Bounds are a box, so use Latin hypercute sampling
    theta_init = KMS_AUX2_drawpoints(multistart_num,dim_p,LB_theta,UB_theta,KMSoptions);
else
    % Bounds define a polytope, so use hit-and-run sampling
    theta_init = KMS_AUX2_drawpoints(multistart_num,dim_p,LB_theta,UB_theta,KMSoptions,A_theta,b_theta,theta_0);
end

% Include theta_0 and bounds
theta_init = [theta_init ; theta_0.'];
M_init = size(theta_init,1);

% Run fminimax over theta_init
theta_fminimax = zeros(M_init,dim_p);
flag_conv      = zeros(M_init,1);
% Objective and constraint:
objective_FeasibleSearch= @(theta)KMS_11_FeasibleSearch_objective(theta,KMSoptions);
constraint_FeasibleSearch = @(theta)KMS_12_FeasibleSearch_constraint(theta,...
            f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,KMSoptions);    

% Run fmincon to find feasible points
if parallel
    parfor ii = 1:M_init
        try
            theta_aug = [theta_init(ii,:).';0];
            [X,~,exitflag] = fmincon(objective_FeasibleSearch,theta_aug,[A_theta, zeros(size(A_theta,1),1)],b_theta,[],[],...
                                [LB_theta;-inf],[UB_theta;inf],constraint_FeasibleSearch,options_fmincon);
            theta_fminimax(ii,:) = X(1:dim_p,1).';
            flag_conv(ii,1) = exitflag;
        catch
            theta_fminimax(ii,:) = theta_init(ii,:);
            flag_conv(ii,1)      = -1;
        end
    end
else
    for ii = 1:M_init
        try
            theta_aug = [theta_init(ii,:).';0];
            [X,~,exitflag] = fmincon(objective_FeasibleSearch,theta_aug,[A_theta, zeros(size(A_theta,1),1)],b_theta,[],[],...
                                [LB_theta;-inf],[UB_theta;inf],constraint_FeasibleSearch,options_fmincon);
            theta_fminimax(ii,:) = X(1:dim_p,1).';
            flag_conv(ii,1) = exitflag;
        catch
            theta_fminimax(ii,:) = theta_init(ii,:);
            flag_conv(ii,1)      = -1;
        end
    end
end

% Keep soltuions that are feasible 
ind = find(flag_conv<= 0);
theta_fminimax(ind,:) = [];

% Check to make sure theta_fminimax(:,ii) is inside the parameter space
A_aug = [A_theta ; eye(dim_p) ; -eye(dim_p)];
b_aug = [b_theta ; UB_theta ; -LB_theta];
size_opt = size(theta_fminimax,1);
ind = find(max(A_aug*(theta_fminimax.') - repmat(b_aug,[1,size_opt])) > 0).';
theta_fminimax(ind,:) = [];

% Compute constraint violation
if ~isempty(theta_fminimax)
    [c_fminimax,CV_fminimax] = KMS_31_Estep(theta_fminimax,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
else
    CV_fminimax =[];
end

% Return of we have found a feasible point:
if ~isempty(find(CV_fminimax == 0))
    flag_feas = 1;
    ind = find(CV_fminimax ==0);
    CV_fminimax = CV_fminimax(ind);
    theta_fminimax = theta_fminimax(ind,:);
    [~,ind_max]  = max(theta_fminimax*p);
    [~,ind_min]  = min(theta_fminimax*p);
    
    theta_feas(1,:) =  theta_fminimax(ind_max,:);
    theta_feas(2,:) =  theta_fminimax(ind_min,:);
    return
end

% If fail to converge, output flag
flag_feas = 0;
theta_feas = [];


end








