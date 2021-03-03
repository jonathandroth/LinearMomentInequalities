function [KMSoptions] = KMSoptions()
%% Code Description: KMS options
% This function outputs a structure called KMSoptions, which are the
% default options for the KMS algorithm.  This structure specifies, among 
% other things, the tuning parameters and maximum iterations in search 
% algorithm.  The user can also add elements to KMSoptions. 
% KMSoptions is passed to, amongother things, the user-specified functions
%
%   moments_w
%   moments_theta
%   moments_gradient
%   moments_stdev
%
% EXAMPLE:
% The user can add the support of covariates to
% KMSoptions:
%   KMSoptions.supp_X = [1 -1 1 -1 ; 1 -1 1 1 ; 1 1 1 -1 ; 1 1 1 1];
% and then call KMSoptions.supp_X in the moments_w function.      

%% User-specified options
% NOTE: If you would like to add to this list, make sure that another
% KMSoption is NOT overwritten.  Otherwise, this program will overwrite it.
KMSoptions.suppX        = [];                                              % Support for X
KMSoptions.psuppX       = [];                                              % Probability for each support point to occur
KMSoptions.DGP          = [];                                              % Specify the DGP
KMSoptions.data_size    = [];                                              % In DGP 2 we need to pass the number of observations to the model-implied moments.
KMSoptions.dX           = [];                                              % Dimension of X (for BCS)
KMSoptions.user_option1 = [];
KMSoptions.user_option2 = [];
KMSoptions.user_option3 = [];
KMSoptions.user_option4 = [];
KMSoptions.user_option5 = [];

%% Key parameters 
KMSoptions.parallel     = 1 ;   % Specify if using parallel computing.  Default is to set parallel = 1 if user has toolbox installed and parallel not specified.
KMSoptions.num_cores    = [];   % Number of cores if using parallel.  If not specified, use max number of cores.
KMSoptions.seed         = 0;    % Seed value  
KMSoptions.CVXGEN       = 1;    % Set equal to 1 if CVXGEN is used.  Set equal to 0 if CVX is used.

% The following parameters can be adjusted by the user in order to improve
% performance/check sensitivity.
KMSoptions.B            = 1001; % Bootstrap repitions                       (Smaller number = faster, but more sensitive.  Recommended set B in [2000,5000].
KMSoptions.EAM_maxit    = 20;   % maximum # of iterations for EAM steps     (Program may terminate early if too small, e.g. less than 20).
KMSoptions.mbase        = 10;   % Initial points in EAM is mbase*dim_p + 1  (Smaller number = faster, but more sensitive.  Recommended values in [10,30])
KMSoptions.h_rate       = 1.8;  % Contraction rate of the parameter space.  (Larger number = faster, but more sensitive.  Recommend values in [1.3,3])
KMSoptions.h_rate2      = 1.25; % Contraction rate for additional points.  Set to a number in (1,h_rate).
KMSoptions.EI_multi_num = 40;   % Maximum number of initial points to use in multi start.  (Smaller = faster, but more sensitive.  Recommend values in [10,100]).
KMSoptions.EAM_obj_tol  = 0.005;% Convergence tolerance.  Once absolute change in objective value < EAM_obj_tol we exit  (provided other condition satisfied)
KMSoptions.EAM_thetadistort = 0.005; % Amount of distrtion of theta#.  We set it equal to EAM_obj_tol in our simulations.

% EAM_maxviol_tol is convergence tolerance level that can be set to either
% inf or a number, such as 1/B.  If EAM_maxviol_tol < infinity, then in
% order to converge we require that the difference between c(theta) and the
% constraint sqrt{n}(f_j + g_j(theta))/sigma_j) is less than
% EAM_maxviol_tol, which we call the maximum violation. 
KMSoptions.EAM_maxviol_tol = inf; 

% How many feasible points to use
% When passing feasible points to EAM, the user can specify to use the
% best avaliable feasible point as a starting point, or can pass all
% avaliable feasible points.
KMSoptions.FeasAll = 0;     % Equal to 0 to pass best, equal to 1 to pass all

% Choose a sampling method.  Set:
% HR = 0 if you want to use uniform and discard method
% HR = 1 if you want to use hit-and-run sampling. 
KMSoptions.HR = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Parameters and options below this line are set to achieve robust
% performance across a wide class of problems.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tolerance and algorithm parameters
% Specify tolerenace and algorithm parameters that
% determine when an optimal or fixed point has been reached.

%%% KMS_0_Main options
% Specify direct_solve = 1 if the user wishes to compute the confidence
% interval using fmincon and solving eq2.16 directly.  
% This method is not recommended as it is slow and may not converge to
% optimum if the moments are not well-behaved.
KMSoptions.direct_solve = 0;

% In the case when we directly solve the problem, we use multistart:
options_multistart = MultiStart;
if KMSoptions.parallel == 1
    options_multistart.UseParallel = 'always';
end
KMSoptions.options_multistart   = options_multistart;
KMSoptions.multistart_num       = 50;
KMSoptions.options_linprog      = optimoptions('linprog','Algorithm','dual-simplex','Display','off'); 

%%% EAM options
% Below are a list of options for the EAM algorithm.
% fmincon options for expected improvement maximization:
KMSoptions.options_fmincon = optimoptions('fmincon','Display','off','GradObj','on','GradConstr','on','Algorithm','sqp','ScaleProblem','none');

% EAM parameters (general/convergence): 
KMSoptions.EAM_tol = 0.0001;                    % Require ||thetaU(iter) -thetaU(iter-1)||<EAM_tol
KMSoptions.EAM_minit = 4;                       % Cannot converge unless we have performed EAM_minit EAM iterations      
% SET ABOVE. KMSoptions.EAM_maxit = 20;         
% SET ABOVE. KMSoptions.mbase = 10;                     
% SET ABOVE. KMSoptions.h_rate = 2;                    
% NB: setting h_rate too large will contract the parameter space too
% quickly, which can bias the confidence interval.  If h_rate is too small,
% then the program might stall. It is recommended setting h_rate between
% 1.4 and 2.  Start with h_rate = 1.4, and if EAM appears fails to update
% the optimal projection for a long string of iterations, increase h_rate
% to 1.6. 

% EAM parameters (E-step): 
% E-step: Fixed-point parameters for finding critical value c
KMSoptions.c_maxit = 100;                           % Maximum number of iterations in root-finding algorithm for critical val
KMSoptions.hc_tol  = 1/KMSoptions.B;                % Tolerence level for root-finding algorithm
KMSoptions.c_tol  = 1e-8;                           % Tolerence level for root-finding algorithm
KMSoptions.Brent_tol  = 1e-5;                       % Tolerence level for root-finding algorithm 
KMSoptions.CVX_maxit = 25;                          % Maximum number of iterations for CVX. 
% NB: CVX usually converges in less than 20 steps and CVXGENs default
% is 25.  We should not set this number too high, because if the set of
% interest is empty, then we run the CVX program equal to the maximum
% number of CVX iterations, and hence consume unneccessary time. 

% EAM parameters (A-step): 
addpath ./dace                                          % use DACE kriging toolbox
KMSoptions.dace_theta       = 10;                       % DACE parameter
KMSoptions.dace_lob         = 1e-1;                     % DACE parameter
KMSoptions.dace_upb         = 20;                       % DACE parameter

% EAM parameters (M-step): 
KMSoptions.unif_points      = 20;                      % Number of uniform points to draw from {theta : p'theta >= p'theta#}
KMSoptions.EI_points        = 10;                      % Number of points with positive expected improvement to draw
KMSoptions.EI_points_start  = 4;                       % Number of points to draw from contracting ball of radius r_min < r < r_max 
KMSoptions.r_min            = 1e-4;                    % Smallest ball to draw positive EI draws from.
KMSoptions.EI_num           = 1;                       % Number of points to keep from expected improvement maximization. 
KMSoptions.unif_num         = 1;                       % Number of points to uniformly draw to add to list on each iteration
% EI_num + unif_num is equal to the number of points to update on each EAM
% iteration, after the first iteration.

% All other parameters (not just EAM)
KMSoptions.siglb            = 1e-4;                     % Lower bound on sigma, the variance parameter                         
KMSoptions.beta             = 0.01;                     % Max bias is equal to beta from searching over rho-box (see Pg 30).
KMSoptions.f_keep_threshold = 1e-4;                     % Threshold for moment functions that are close to boundary.

%%  Placeholders -- DO NOT CHANGE
KMSoptions.p          = [];      
KMSoptions.dim_p 	  = [];
KMSoptions.n          = [];
KMSoptions.rho        = [];
KMSoptions.J          = [];
KMSoptions.J1         = [];
KMSoptions.J2         = [];
KMSoptions.J3         = [];
KMSoptions.S          = [];
KMSoptions.paired_mom = [];
KMSoptions.c_B        = [];
KMSoptions.kappa      = [];
KMSoptions.phi        = [];
KMSoptions.alpha      = [];
KMSoptions.LB_theta   = [];
KMSoptions.UB_theta   = [];
KMSoptions.A_theta    = [];
KMSoptions.b_theta    = [];
KMSoptions.CI_lo      = [];
KMSoptions.CI_hi      = [];
KMSoptions.theta_lo      = [];
KMSoptions.theta_hi      = [];
KMSoptions.sample_method = [];                                                  
KMSoptions.type          = [];
KMSoptions.e_points_init = [];
KMSoptions.CVXGEN_name   = [];
KMSoptions.CI_method     = [];
KMSoptions.component    = [];


end
