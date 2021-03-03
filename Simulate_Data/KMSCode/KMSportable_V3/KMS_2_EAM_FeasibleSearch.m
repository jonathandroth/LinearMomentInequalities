function [theta_feas,flag_feas,theta_out,c_out,CV_out,maxviol_out] =  KMS_2_EAM_FeasibleSearch(p,theta_0,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions)
%% Code description: EAM Feasible Search
%  This function executes an auxiliary search for a feasible point.
%  The algorithm attempts to find a feasible point by minimizing the 
%  "hard moment" constraint:
%
%    min max_{j=1,...,J} (sqrt(n)m(X,theta)/sigma(X) - c_L(theta))
%
% where c_L(theta) is approxiamted via DACE.  
%
% In the case that the optimal point from this search satisfies 
%
%    CV = sum_j max(0, (sqrt(n)m(X,theta)/sigma(X) - c(theta))^2) = 0
% 
% we take theta_feas to be this feasible point. If the point that solves
% this problem is not feasible, then we add it to the list of evaluation
% points and re-update the DACE model for c_L(theta), as well as a uniformly
% drawn point.  We run this procedure until either we find a feasible point
% or we reach the maximum number of iterations.  If, after this, we fail to
% find a point, we conclude that no feasible point exists.  This could be
% due to model miss-specificaiton, a bug in the user-specified function,
% or the model does not satisfy regularity conditions. Alternatively,  
% the feasible search method may not work due to instability in the
% gradient function  (or another function).
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
e_points_init      = KMSoptions.e_points_init;
dace_theta         = KMSoptions.dace_theta;
dace_lob           = KMSoptions.dace_lob;
dace_upb           = KMSoptions.dace_upb;
EAM_maxit          = KMSoptions.EAM_maxit;

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

% Draw points from Theta for initial DACE approx.
if sample_method == 0 || sample_method == 1 
    % Bounds are a box, so use Latin hypercute sampling
    theta_init = KMS_AUX2_drawpoints(e_points_init,dim_p,LB_theta,UB_theta,KMSoptions);
else
    % Bounds define a polytope, so use hit-and-run sampling
    theta_init = KMS_AUX2_drawpoints(e_points_init,dim_p,LB_theta,UB_theta,KMSoptions,A_theta,b_theta,theta_0);
end

% Set theta_Estep (initial set of points for EAM) to be the union of the
% uniformly drawn points, bounds on Theta, and theta_feas.  
theta_Estep = [theta_init ; theta_0.' ];

% theta_Astep is the set of thetas that EAM has explored
% So that the program runs correctly, I initiate it to be the empty set.
% theta_Astep will be updated shortly in the EAM algorithm
theta_Astep = [];
c_Astep = [];
CV_Astep =[];
maxviol_Astep = [];

%% EAM Feasible Search Routine
for iter=1:EAM_maxit
    % Step 1) E-step
    [c_Estep,CV_Estep,maxviol_Estep] = KMS_31_Estep(theta_Estep,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
   
    % Update (theta,c) for the A-step
    theta_Astep = [theta_Astep ; theta_Estep];
    c_Astep     = [c_Astep ; c_Estep];
    CV_Astep = [CV_Astep;CV_Estep];
    maxviol_Astep = [maxviol_Astep;maxviol_Estep];
    
    % Keep only unique points
    [theta_Astep,ind] =  unique(theta_Astep,'rows'); 
    c_Astep = c_Astep(ind);
    CV_Astep = CV_Astep(ind);
    maxviol_Astep = maxviol_Astep(ind);
    
    % Step 1.5) Check feasibility
    feas = find(CV_Astep == 0);
    if ~isempty(feas) 
        theta_temp = theta_Astep(feas,:);
        [~,ind_max] = max(theta_temp*p);
        [~,ind_min] = min(theta_temp*p);
        theta_feas = theta_temp([ind_max;ind_min],:).';
        theta_feas  =  unique(theta_temp,'rows'); 
        theta_out   = theta_Astep;
        c_out       = c_Astep;
        CV_out      = CV_Astep;
        maxviol_out = maxviol_Astep;
        flag_feas  = 1;
        return
    end
    
    % Step 2) A-Step 
    % This step interpolates critical values outside of the evaluation 
    % points
    [theta_dmodel,ind] = uniquetol(theta_Astep,1e-10,'ByRows',true); 
    c_dmodel = c_Astep(ind,:);
    dmodel = dacefit(theta_dmodel,c_dmodel,@regpoly0,@corrgauss,dace_theta,dace_lob,dace_upb);
    
    % Step 3) M-Step 
    % This step draws new point for next iteration
    
    % Draw initial points for fminimax
    % Draw points from Theta and in a parfor loop pass these to fmincon.
    if sample_method == 0 || sample_method == 1 
        % Bounds are a box, so use Latin hypercute sampling
        theta_fminimax = KMS_AUX2_drawpoints(multistart_num,dim_p,LB_theta,UB_theta,KMSoptions);
    else
        % Bounds define a polytope, so use hit-and-run sampling
        theta_fminimax = KMS_AUX2_drawpoints(multistart_num,dim_p,LB_theta,UB_theta,KMSoptions,A_theta,b_theta);
    end
    
    % Include theta_0 and bounds
    theta_fminimax = [theta_fminimax ; theta_0.' ];
    num_Mstep = size(theta_fminimax,1);
    
    % Run fminimax over theta_init
    theta_Mstep = zeros(num_Mstep,dim_p);
    gamma_Mstep = zeros(num_Mstep,1);
    flag_conv   = zeros(num_Mstep,1);
    
    % Objective and constraint:
    objective_FeasibleSearch= @(theta)KMS_21_EAM_FeasibleSearch_objective(theta,KMSoptions);
    constraint_FeasibleSearch = @(theta)KMS_22_EAM_FeasibleSearch_constraint(theta,...
            f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,dmodel,KMSoptions);   
    
    % Find solution(s)
    if parallel
        parfor ii = 1:num_Mstep
            try
                theta_aug = [theta_fminimax(ii,:).';0];
                [X,fval,exitflag] = fmincon(objective_FeasibleSearch,theta_aug,[A_theta, zeros(size(A_theta,1),1)],b_theta,[],[],...
                                    [LB_theta;-inf],[UB_theta;inf],constraint_FeasibleSearch,options_fmincon);
                theta_Mstep(ii,:) = X(1:dim_p,1).';
                gamma_Mstep(ii,1) = fval;
                flag_conv(ii,1)     = exitflag;
            catch
                theta_Mstep(ii,:)   = theta_fminimax(ii,:);
                gamma_Mstep(ii,1)   = 0;
                flag_conv(ii,1)     = -1;
            end
        end
    else
        for ii = 1:num_Mstep
            try
                theta_aug = [theta_fminimax(ii,:).';0];
                [X,fval,exitflag] = fmincon(objective_FeasibleSearch,theta_aug,[A_theta, zeros(size(A_theta,1),1)],b_theta,[],[],...
                                    [LB_theta;-inf],[UB_theta;inf],constraint_FeasibleSearch,options_fmincon);
                theta_Mstep(ii,:) = X(1:dim_p,1).';
                gamma_Mstep(ii,1) = fval;
                flag_conv(ii,1)     = exitflag;
            catch
                theta_Mstep(ii,:)   = theta_fminimax(ii,:);
                gamma_Mstep(ii,1)   = 0;
                flag_conv(ii,1)     = -1;
            end
        end
    end 
     % Keep solutions that are feasible 
    ind = find(flag_conv<= 0);
    theta_Mstep(ind,:) = [];
    gamma_Mstep(ind,:) = [];
    
    % Check solutions are in the parameter space
    A_aug = [A_theta ; eye(dim_p) ; -eye(dim_p)];
    b_aug = [b_theta ; UB_theta ; -LB_theta];
    size_opt = size(theta_Mstep,1);
    ind = find(max(A_aug*(theta_Mstep.') - repmat(b_aug,[1,size_opt])) > 0).';
    theta_Mstep(ind,:) = [];
    gamma_Mstep(ind,:) = [];

    % Sort by gamma_Mstep
    [gamma_Mstep,I] = sort(gamma_Mstep,'ascend');
    theta_Mstep = theta_Mstep(I,:);
    
    % Update theta_Estep by picking minimizer + a uniform point
    if isempty(I) == 0
        theta_Estep = theta_Mstep(1,:);
    else
        theta_Estep = [];
    end
    
    if sample_method == 0
       theta_draw = KMS_AUX2_drawpoints(1,dim_p,LB_theta,UB_theta,KMSoptions);
    else
        theta_draw = KMS_AUX2_drawpoints(1,dim_p,LB_theta,UB_theta,KMSoptions,A_theta,b_theta);
    end
    theta_Estep = [theta_Estep ; theta_draw];
    
end

% If failed to converge, output empty and flag
theta_feas = [];
flag_feas = 0;
theta_out   = [];
c_out       = [];
CV_out      = [];
maxviol_out = [];
end








