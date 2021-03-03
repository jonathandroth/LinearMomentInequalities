function [KMS_confidence_interval,KMS_output] = KMS_0_Main(W,theta_0,p,theta_feas,LB_theta,UB_theta,A_theta,b_theta,alpha,type,CI_method,kappa,phi,CVXGEN_name,KMSoptions)
%% Code Description: Main File
% CODE AUTHORS:
% 	Hiroaki Kaido
%   Francesca Molinari
% 	Joerg Stoye
% 	Matthew Thirkettle  
% 
% REFERENCE PAPER:
% 	Confidence Intervals for Projections of Partially Identified Parameters
% 	By Hiroaki Kaido, Francesca Molinari, and Joerg Stoye
%   CeMMap Working Paper CWP02/16.
% 
% VERSION | DATE       | MODIFICATION
% -------------------------------------------------------------------------
% V*      | 10-03-2017 | (Briefly describe modification, being specific)
%         |            |
% -------------------------------------------------------------------------
% V1      |            | p can be any directional vector
%         |            | Included polytope constraints A_theta*theta <=b_theta
%         |            | Included an additional feasible search if the
%         |            | first fails that is based on EAM.
% -------------------------------------------------------------------------
% V2      | 22-08-2017 | BUG FIXES:
%         |            | (KMS_3_EAM, 215-219) If parameter is too close to
%         |            |  the bounday, due to numerical error the LB or UB
%         |            |  can be overridden with a value outside of the
%         |            |  user-specified LB/UB causing an error in a
%         |            |  function.  
%         |            | (KMS_0_Main, 477-479) The parallel bootstrap caused a random
%         |            |  subseed to be generated.  The seed is fixed after
%         |            |  bootstrap.  
%         |            |  (KMS_3_EAM, 312-345) In DGP = 7, fmincon
%         |            |  attempted to evaluate a point outside of the 
%         |            | user-specified bounds.  For now we put a try/catch
%         |            | command in.  Will investigate later.
%         |            | (KMS_AUX2_drawpoints, 57-99) Failed to evaluate S
%         |            | if S is empty due to previous failure.  Update
%         |            | S iff hit and run did not fail.  When hit-and-run
%         |            | fails, often lcon2vert fails as well.  This is because
%         |            | the parameter is close to the boundary and the
%         |            | effective search-set is thin. Both hit-and-run and 
%         |            | lcon2vert require a non-thin set.  For now we
%         |            | a try/catch and return S = [];.  If S = [], then
%         |            | we likely only use a uniform draw.
%         |            | (lcon2vert) Replaced keyboards with error.
% -------------------------------------------------------------------------
% V3      | 12-10-2017 | Due to numerical problems with hit-and-run sampling, 
%         |            | we implemented the draw-and-discard method for
%         |            | DGP8.  Method currently only works if p is a bais
%         |            | vector.
% -------------------------------------------------------------------------
% 
% SUMMARY:
% This program executes the KMS procedure to find the confidence interval
% of a parameter theta that is identified by moment inequalities and moment
% equalities:
% 
%       E_P[m_i(W,theta)] <= 0 for i = 1,...,J1
%       E_P[m_i(W,theta)] = 0  for i = J1+1,...,J1+J2
% 
% The output is an interval [thetaL,thetaU] that covers a
% particular direction p of theta with (1-alpha) probability.
% 
% This particul function assumes that all moment (in)equalities are 
% separable in data W and parameter theta.  That is, there exists functions
% f,g such that
% 
%   E_P[m_i(W,theta)]  =E_P[f_i(W)] + g_i(theta) for all i = 1,...,J, J =J1+J2
% 
% Code for the non-separable case will be included in a future release.
% For the best performance, it is recommended that the user writes moments 
% in separable form if possible.  
% 
% This code has also been optimized for parallel computing.
% 
% CVXGEN INSTRUCTIONS
% We use the high-speed quadratic solver CVXGEN to find whether or not the
% constraint set intersects the hyperplane defined by p'*lambda
% (See Pg11-12, eq 2.11-2.12).  The solver CVX can also be used instead.
% CVXGEN is appropriate when the LP problem is small. For each problem a
% new CVXGEN MEX file needs to be compiled.  See manual for instructions.
% 
% USER-SPECIFIED FUNCTIONS
% Among other things, the user specifies the moment conditions,
% gradients of the moments, and standard deviation estimator of moments.
% The first step is to write the problem in separable form:
% 
%   E_P[m(W_i,theta)]  = E_P[f(W_i)] + g_i(theta).
% 
% Next, specify the following functions:
%       - f(W_i).      This is done by modifying moments_w.
%       - g(theta).    This is done by modifying moments_theta.
%       - Dg(theta).   This is done by modifying moments_gradient.
%       - sigma(W_i).    This is done by modifying moments_stdev.
% Further instructions on this is provided in each of these functions and
% in the manual.
% 
% Last, specify input variables for this function (KMS_0_Main)
% 
% INPUTS:
%   Data:
%   (1) W.          Data vector W = [w_1, w_2, ... , w_K], where w_k is n-by-1.  
%                   n is the sample size.
% 
%   Inputs for the parameter theta:
%   (2) theta_0.    This is initial guess for theta.
%                   theta_0 is dimension dim_p-by-1.
% 
%   (3) p.          This is the directional vector that is dim_p-by-1.  
%                   We solve min/max p'theta s.t. normalized moments less 
%                   than or equal to c(theta).  In order to get a confidence
%                   set for a particular parameter,say theta_k, set p to be
%                   the basis vector.
% 
%   (4) theta_feas. K-by-dim_p matrix of feasible points. We perform an 
%                   auxiliary for points that satisfy 
%                   bar m_j(theta)]/sigma_j <= c(theta)  for all j. 
%                   We provide the option for the user to provide points
%                   that satisfy this condition. See user-notes below.
% 
%   (5)-(8) LB_theta,UB_theta,A_theta,b_theta. 
%                   This defines the parameter space.  It is required that 
%                   the parameter space is bounded.  The parameter space is
%        {theta : LB_theta <= theta <= UB_theta & A_theta*theta <= b_theta}
%                   It is required that LB_theta and UB_theta are specified.
%                   The polytope constraints A_theta and b_theta can be set
%                   to the emptyset, in which case the parameter space is
%                   the simple box. LB_theta and UB_theta are dim_p-by-1,
%                   A_theta is L-by-dim_p (L determined by user), and 
%                   b_theta is L-by-1.
% 
%   (9) alpha.      This is the nominal level of converage.  Valid input is
%                   alpha in (0,1).  Recommended levels include 
%                   alpha = 0.1, 0.05 and 0.01 for the 90%, 95%, and 99% 
%                   confidence intervals.
% 
%   (10) type.      type = 'two-sided' or 'one-sided-LB' or  'one-sided-UB'
%                   to compute a two-sided or one-sided confidence interval
%                   NB: one-sided-LB outputs [-s(-p,C),UB) and one-sided-UB 
%                   outputs (LB,s(p,C)], where UB and LB are the upper and 
%                   lower bounds on the parameter space in direciton p.
% 
%   (11) method.    Equal to 'KMS' or 'AS', to specify the calibrated
%                   or AS projection interval
% 
%   Tuning parameters and procedural functions
%   (12)            kappa.  This is the GMS tuning parameter.  Set
%                   kappa = @(n)kappa_function(n) if user chooses to supply
%                   own function.  Set kappa = NaN to use default parameter
%                   sqrt(log(n)) (see Pg 9, eq 2.8).
% 
%   (13) phi.       This is the GMS function.  Set
%                   phi = @(x)phi_function if user chooses to supply own
%                   function. Set phi = NaN for hard threshold GMS.
%                   (See Pg 9, eq 2.9.)
% 
%   (14) CVXGEN_name. 
%                   This is the name of the csolve.mexw64 name.  If not
%                   specified, will use 'csolve.mexw64'. 
% 
%   (15) KMSoptions.   
%                   This is a structure of additional inputs held
%                   constant over the program.  In the 2x2 entry game, 
%                   KMSoptions includes the support for the covariates and
%                   the probability of support point occuring.  There are 
%                   also options in KMSoptions to  specify  optimization 
%                   algorithm, tolerance, and tuning parameters. 
% 
% 
% % NOTES TO USER:
% 1) The hit-and-run method is sometimes unreliable.  Caution is 
% recommended if polytope constraints are set (A_theta,b_theta) or if p is
% not a basis vector
% 2) In order to reliably implement the EAM algorithm, we require a
% feasible parameter -- a theta that satisfies the moment constraints less than
% or equal to the critical value c(theta).  We implement a pre-search for a
% feasible point.  The ease of finding a starting point depends on the
% specific application, and a larger number of starting draws may be needed
% if there are a lot of inequalities and if the identified set is
% relatively small.  In addition, we allow the user to pass a feasible
% point (found via a user-implemented pre-search).  In particular, if
% theta_feas (an input to this function) is not specified and left empty,
% then our presearch is executed.   Otherwise, the presearch is bypassed
% and the user-specified theta_feas is used instead.
% 
% Feasibility is determined by the function KMS_31_Estep.  The evaluation
% point theta_Estep is feasible iff CV_Estep is equal to 0.  
% Therefore, if theta_feas is provided, it is required that that
% CV_Estep = 0. An error occurs otherwise.

%% Extract relevant information from KMSoptions
B                   = KMSoptions.B;                 % Bootstrap repetitions
beta                = KMSoptions.beta;              % bias parameter
parallel            = KMSoptions.parallel;          % =1 iff parallel computing set
num_cores           = KMSoptions.num_cores;         % cores to use
options_fmincon     = KMSoptions.options_fmincon;   % fmincon options
options_multistart  = KMSoptions.options_multistart;% Multistart options
seed                = KMSoptions.seed;              % Seed number for reproducibility
siglb               = KMSoptions.siglb;
multistart_num      = KMSoptions.multistart_num;

%% Model parameters
% Compute model parameters implied by the user-specified input parameters.
dim_p = size(theta_0,1);     % Dimension of parameter of interest
n     = size(W,1);           % Number of observations
if size(p) ~=  size(theta_0)
    error('Invalid input for p.')
end
% Normalize p
p  = p/norm(p,2);

% Empirical moments:
% Estimator for E_P[f(W_i)], where E_P[m(W_i,theta)] = E_P[f(W_i)] + g(theta).
% Emprical moments are obtained first, as this defines the number of moment
% conditions
[f_ineq,f_eq,f_ineq_keep,f_eq_keep,paired_mom,J1,J2,J3]  = moments_w(W,KMSoptions);

% Total number of moments:
J =  J1 + 2*J2;                                                   
if J1 < 0 || J2 < 0 || J3 < 0 || J <= 0
    error('Invalid value of J.')
end
rho = -norminv((1-(1-beta).^(1/(dim_p*nchoosek(J,dim_p))))/2);      % size of rho-box (see Pg 30).
c_B = norminv(1-alpha/J);                                           % Bonferroni critical value (Pg 10, after Eq 2.11)

%% Error checks and adjustments
% Check to make sure inputs are valid
if size(theta_0,2) > 1
    error('Invalid dimension for theta.  It should be dim_p-by-1.')
end

if strcmp('one-sided-LB',type) == 0 &&  strcmp('one-sided-UB',type) == 0 && strcmp('two-sided',type) ==0
    error('Invalid input for type.')
end

if ~isempty(theta_feas)
    if size(theta_feas,2) ~= dim_p
         error('Invalid input for theta_feas.  Set theta_feas to emptyset or K-by-dim_p matrix of feasible points')
    end
end

try isnan(kappa);
    kappa = @(x)sqrt(log(x));
catch
    try kappa(n);    
    catch
        error('Input for kappa not valid.  Either set kappa = nan or kappa to a function handle with one input.')
    end
end
if kappa(n) <= 0
     error('Input for kappa not valid.  kappa function needs to satisfy kappa(n) > 0 and kappa(n) -> infinity.')
end

% GMS function
if isnan(phi)
    phi = @(x)KMS_AUX1_phifunc(x);
elseif ~(strcmp(class(phi),'function_handle') && nargin(phi) == 1)
     error('Input for phi not valid.  Either set phi = 1 for hard thresholding or phi to a function handle with one input.')
end

% Parameter space
if isempty(LB_theta) || isempty(UB_theta)
    error('Box constraints on the parameter space are required.')
end
if ~isempty(LB_theta) && max(size(LB_theta) ~= [dim_p,1])
   error('Lower bound on parameter space incorrect size.  LB_theta should be either dimension dim_p-by-1 or empty.')
end
if ~isempty(UB_theta) && max(size(UB_theta) ~= [dim_p,1])
    error('Upper bound on parameter space incorrect size. UB_theta should be either dimension dim_p-by-1 or empty.')
end

% Number of linear constraints on parameter space:
L = size(A_theta,1); 
if ~isempty(A_theta) && size(A_theta,2) ~= dim_p
    error('Linear constraints on parameter space incorrect size. A_theta should be dimension L-by-dim_p.')
end
if ~isempty(b_theta) && max(size(b_theta) ~= [L,1])
    error('Linear constraints on parameter space incorrect size. b_theta should be dimension L-by-1.')
end
if (isempty(A_theta) && ~isempty(b_theta)) ||  (~isempty(A_theta) && isempty(b_theta))
    error('Require input of either both A_theta and b_theta or neither.')
end

% Determine sampling method
% We initially sample from the parameter space uniformly in order to obtain
% initial points for the EAM algorithm.  If the parameter space is a
% polytope, we use hit-and-run sampling or draw-and-discard sampling. 
% Otherwise, we use uniform sampling.
% We also contract the parameter space.  The contraction direction is
% determined by the directional vector p.  If p is a basis vector and the
% parameter space is determined by box constraints, we can use simple
% uniform sampling.  Otherwise, we use HR or DD sampling.
if isempty(A_theta) && abs(max(p)) == 1
    % p is a basis vector and the parameter space is a simple box.  So
    % always use uniform sampling.
    sample_method = 0;
    [~,component] = max(abs(p));
    KMSoptions.component = component;
elseif isempty(A_theta)
    % p is not a basis vector, but the parameter space is a simple box.
    % Use uniform sampling initially, then switch to HR sampling 
    % in the contraction phase.
    sample_method = 1;
elseif KMSoptions.HR == 0
    % Use DD method
    sample_method = 3;
else
    % Use HR method
    sample_method = 2;
end

% Augment constraints
A_aug = [A_theta ; eye(size(UB_theta,1)) ; -eye(size(LB_theta,1))];
b_aug = [b_theta ; UB_theta   ; -LB_theta];
L_aug   = size(b_theta,1);

% Check to see if theta_0 is feasible
if min(A_aug*theta_0 <= b_aug) == 0
    error('theta_0 is not feasible.')
end

% Check to see if the parameter space is bounded in all directions
options_linprog = KMSoptions.options_linprog;
for ii = 1:dim_p
    direction = zeros(dim_p,1);
    direction(ii,1) = 1;
    [~,~,flag_space] =  linprog(direction,A_aug,b_aug,[],[],[],[],[],options_linprog);
    if flag_space ~= 1
         error('Parameter space is unbounded: re-adjust either the upper/lower bound or linear constraints A_theta,b_theta.')
    end
    direction(ii,1) = -1;
    [~,~,flag_space] =  linprog(direction,A_aug,b_aug,[],[],[],[],[],options_linprog);
    if flag_space ~= 1
         error('Parameter space is unbounded: re-adjust either the upper/lower bound or linear constraints A_theta,b_theta.')
    end
end

% Determine the min/max value of p'theta s.t. theta is in the parameter
% space
[theta_lo,CI_lo,~] =  linprog(p,A_aug,b_aug,[],[],[],[],[],options_linprog);
[theta_hi,CI_hi,~] =  linprog(-p,A_aug,b_aug,[],[],[],[],[],options_linprog);
CI_hi = -CI_hi;

% Pre-specify number of cores to be used.  Truncate if this number exceeds
% number of cores on the computer/server.  Set num_cores = 1 if parellel
% toolbox is not installed.
v_=ver;
[installedToolboxes{1:length(v_)}] = deal(v_.Name);
parallel_installed = all(ismember('Parallel Computing Toolbox',installedToolboxes));
parallel  = min(parallel, parallel_installed);
KMSoptions.parallel = parallel;
if parallel == 0
    % Do nothing
elseif isempty(num_cores)
    % If num_cores not specified, set num_cores to max
    num_cores = feature('numCores');
    KMSoptions.num_cores = num_cores;
else
    % Otherwise, set num_cores to min of specified cores and maximum
    % possible cores
    num_cores = min(num_cores,feature('numCores'));
    KMSoptions.num_cores = num_cores;
end

% Set CVX_Name to default if not specified
if isempty(isempty(CVXGEN_name))
    CVXGEN_name = 'csolve';
end

% Boot up parallel toolbox
if isempty(gcp('nocreate')) && parallel ==1
    % Create parallel pool
    parpool(num_cores);
elseif parallel ==1
    % Otherwise, check to make sure correct number of workers booted
    % If not, delete pool and create a new one with correct number of 
    % workers.  p = gcp
    par_pool_gcp = gcp;
    if par_pool_gcp.NumWorkers ~= num_cores
        delete(gcp)
        parpool(num_cores);
    end
end

%% Additional parameters to add to KMSoptions
addpath ./CVXGEN                                                % MEX files for CVXGEN
addpath ./dace                                                  % DACE parameter (A-step in EAM)
addpath ./hitandrun                                             % Sampling method
addpath ./polytopes/                                            % Vertice sampling
addpath ./MVNorm                                                % Multivariate norm c code 
KMSoptions.dace_theta = KMSoptions.dace_theta*ones(1,dim_p);    % DACE parameter (A-step in EAM)
KMSoptions.dace_lob   = KMSoptions.dace_lob*ones(1,dim_p);      % DACE parameter (A-step in EAM)
KMSoptions.dace_upb   = KMSoptions.dace_upb*ones(1,dim_p);      % DACE parameter (A-step in EAM)
KMSoptions.p          = p;                                      % Projection direction
KMSoptions.dim_p 	  = dim_p;                                  % Dimension of parameter space
KMSoptions.n          = n;                                      % Sample size
KMSoptions.rho        = rho;                                    % Size of rho-box
KMSoptions.J          = J;                                      % Total number of moments
KMSoptions.J1         = J1;                                     % Number of moment inequalities
KMSoptions.J2         = J2;                                     % Number of moment equalities
KMSoptions.J3         = J3;                                     % Number of paired moment inequalities
KMSoptions.paired_mom = paired_mom;                             % Vector identifying paired moment inequalities
KMSoptions.c_B        = c_B;                                    % First-pass upper bound on critical value c(theta)
KMSoptions.kappa      = kappa(n);                               % GMS tuning parameter
KMSoptions.phi        = phi;                                    % GMS function
KMSoptions.alpha      = alpha;                                  % Significance level
KMSoptions.LB_theta   = LB_theta;                               % Parameter space
KMSoptions.UB_theta   = UB_theta;                               % Parameter space
KMSoptions.A_theta    = A_theta;                                % Linear constraint on parameter space
KMSoptions.b_theta    = b_theta;                                % Linear constraint on parameter space
KMSoptions.CI_lo      = CI_lo;                                  % Lower bound on projected confidence interval implied by parameter space
KMSoptions.CI_hi      = CI_hi;                                  % Upper bound on projected confidence interval implied by parameter space
KMSoptions.theta_lo   = theta_lo;                               % Corresponding vector to CI_lo
KMSoptions.theta_hi   = theta_hi;                               % Corresponding vector to CI_hi
KMSoptions.sample_method = sample_method;                       % Indicates which sampling method to use.
KMSoptions.type          = type;                                % One-sided or two-sided test
KMSoptions.CVXGEN_name   = CVXGEN_name;                         % Name of user-built CVXGEN .mex file
KMSoptions.e_points_init = dim_p*KMSoptions.mbase + 1;          % # of initial evaluation points
KMSoptions.CI_method     = CI_method;                           % Confidence interval method: KMS or AS 


%% Re-centered Empirical Process
% (See Pg 8, eq 2.6)
% We compute the re-centered empirical process G.  This can be done
% outside of the EAM loop in the separable case.  The bth recentered
% empirical process is:
%
%       G_{ib} = sqrt(n)( (1/n) sum_i (f_j(W_i^b) -  f_j)/sigma_j
%
% for b=1,...,B.  G_{ib} is moment i evaluated at the bootstrapped data W^b
% subtract the the empirical mean of the bootstrapped data,
% (1/n)sum_if_j(W_i), normalized by the standard deviation.
%
% Note that the mean of the empirical distribution has already been
% obtained, and is equal to f_ineq, f_eq.
%
% BOOTSTRAP:
% Input is the data matrix W.  We output moments
% computed on bootstrapped data
if J1 > 0 && J2 > 0
    f_ineq_b = zeros(J1,B);
    f_eq_b = zeros(2*J2,B);
elseif J1 > 0 && J2 == 0
    f_ineq_b = zeros(J1,B);
    f_eq_b = [];
else
    f_ineq_b =[];
    f_eq_b = zeros(2*J2,B);
end
if parallel == 1
    parfor bb = 1:B
        % Set subseed so that we can reproduce parallelized output
        stream=RandStream('mlfg6331_64','Seed',seed);
        RandStream.setGlobalStream(stream);
        stream.Substream = bb + B*10^3;

        % Resample W
        % First, draw an set of indices of length n (= number of
        % observations).  This will be used to select observations (rows) from
        % W for our bootstrap sample.
        ind = ceil(n*rand(n,1));
        W_b = W(ind,:);
        if J1 > 0 && J2 > 0
            [f_ineq_b(:, bb),f_eq_b(:, bb)] = moments_w(W_b,KMSoptions);
        elseif J1 > 0 && J2 == 0
             [f_ineq_b(:, bb),~] = moments_w(W_b,KMSoptions);
        else
            [~,f_eq_b(:, bb)] = moments_w(W_b,KMSoptions);
        end
    end
else
    for bb = 1:B
        % Set subseed so that we can reproduce parallelized output
        stream=RandStream('mlfg6331_64','Seed',seed);
        RandStream.setGlobalStream(stream);
        stream.Substream = bb + B*10^3;

        % Resample W
        % First, draw an set of indices of length n (= number of
        % observations).  This will be used to select observations (rows) from
        % W for our bootstrap sample.
        ind = ceil(n*rand(n,1));
        W_b = W(ind,:);
        if J1 > 0 && J2 > 0
            [f_ineq_b(:, bb),f_eq_b(:, bb)] = moments_w(W_b,KMSoptions);
        elseif J1 > 0 && J2 == 0
             [f_ineq_b(:, bb),~] = moments_w(W_b,KMSoptions);
        else
            [~,f_eq_b(:, bb)] = moments_w(W_b,KMSoptions);
        end
    end
end

% Fix substream -- if ran in parallel, substream is incorrect
stream=RandStream('mlfg6331_64','Seed',seed);
RandStream.setGlobalStream(stream);
stream.Substream = B + B*10^3 + 1;

% STANDARD DEVIATION
% Compute the estimate for the standard deviation 
[f_stdev_ineq,f_stdev_eq]= moments_stdev(W,f_ineq,f_eq,J1,J2,KMSoptions);

% Truncate standard errrors
f_stdev_ineq = max(siglb,f_stdev_ineq);
f_stdev_eq = max(siglb,f_stdev_eq);

% RECENTER BOOTSTRAP MOMENTS
G_ineq = sqrt(n).*(f_ineq_b - repelem(f_ineq,1,B))./repelem(f_stdev_ineq,1,B);
G_eq   = sqrt(n).*(f_eq_b - repelem(f_eq,1,B))./repelem(f_stdev_eq,1,B);

%% Feasibility
% If the user specifies a set of feasible points, theta_feas, we check that
% these points are feasible for the problem.  We also order the
% points from maximizer to minimizer of p'theta.
if ~isempty(theta_feas)
    [~,CV_feas,~] = KMS_31_Estep(theta_feas,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
    theta_feas = theta_feas(CV_feas==0,:);
    if isempty(theta_feas)
        error('User provided matrix of feasible points, theta_feas, none of which are feasible.')
    end
    
    % Order theta_feas from largest value of p'theta to lowest value of p'theta
    val_feas = theta_feas*p;
    [~,I] = sort(val_feas,'descend');
    theta_feas = theta_feas(I,:);
end

%% EAM - Jones' Method
% The EAM algorithm is based on Jones' Method.  The algorithm is executed
% as follows:
%
% The output is the optimal upper and lower bound that max/min dot(p,theta)
% subject to being in the constraint set.  On the first iteration we:
%
%   (Initialize)  Draw uniformly theta from the parameter space Theta,
%      say theta_1,...,theta_K.  Also perform an initial search that
%      minimizes the standardized moments.  This is done in order to find
%      at least one feasible point with respect to Eq 2.16 on Pg 13.
%   Four steps per iteration:
%   1) Evaluate:
%        Compute the critical value c(theta) at each theta_l, as well as
%        the constraint violation WRT eq 2.16 on Pg 13.
%   2) Approximate:
%        Use auxiliary model to interpolate c(theta)
%   3) Maximize:
%        Find theta* that maximizes expected improvement (EI).  The EI
%        maximizer is a point that improves the objective function p'theta
%        that is likely (determined by a Gaussian prior on the constraint)
%       to satisfy the constraint 
%       normalized moments <= critical value c(theta)
%   4) Convergent criteria:
%   Check to see if stopping conditions are satisfied, which requires that
%   the step-size of the optimal value is small and expected improvement
%   small relative to the optimal value dot(p,theta).  
%
% We run the EAM algorithm in direction p and -p to find the upper and
% lower bounds, respectively.

% Auxilary Search For Feasible Point
% We first perform an auxilary search to find a feasible point.  The same
% point can be used for the upper and lower bound.  If we find more than
% one feasible point we will select the best point according to the rule:
%   min/max p'theta
if isempty(theta_feas)
    [theta_feas,flag_feas]  = KMS_1_FeasibleSearch(p,theta_0,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
else
    flag_feas = 1;
end
% In the case when flag_feas = 0, we try an EAM-style search for a feasible
% point.  To do this, we first approximate c(theta) using c_L(theta).  We
% then execute a feasible search similar to KMS_1_EAM_FeasibleSearch,
% except with constraints m_j(X,theta) <= c_L(theta).  On each iteration,
% we check to see if the solution(s) is feasible.  If not feasible, we add
% it to the set of evaluation points.
if flag_feas == 0
    [theta_feas,flag_feas,theta_init,c_init,CV_init,maxviol_init]  = KMS_2_EAM_FeasibleSearch(p,theta_0,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
else
    theta_init = [];
    c_init     = [];
    CV_init    = [];
    maxviol_init = [];
end

% If we still cannot find a solution, we output an empty set for the KMS
% interval.  If we cannot find a solution, either the model is 
% miss-specified or there is a bug in the user-written command for
% m(X,theta) or Dm(X,theta).  
if flag_feas == 0
    warning('Cannot find feasible point -- Output empty set for KMS interval')
    KMS_output.flag_feas = 0;
    KMS_confidence_interval = [CI_lo,CI_hi];
    return
else
    KMS_output.flag_feas = 1;
end

% Save feasible to KMS_output
KMS_output.thetafeasU = theta_feas(1,:);
KMS_output.thetafeasL = theta_feas(1,:);

% EAM Algorithm
if ((strcmp('one-sided-UB',type) == 1 || strcmp('two-sided',type) ==1)) && flag_feas == 1
    % Find maximum
    tU_EAM = tic;
    if KMSoptions.FeasAll
        [thetaU_hat,thetaU_optbound,c,CV,EI,flagU_opt]  = KMS_3_EAM(p,1,theta_feas,theta_init,c_init,CV_init,maxviol_init,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
    else 
        [thetaU_hat,thetaU_optbound,c,CV,EI,flagU_opt]  = KMS_3_EAM(p,1,theta_feas(1,:),theta_init,c_init,CV_init,maxviol_init,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
    end
    % Save additional output to the KMS_output structure
    KMS_output.thetaU_EAM = thetaU_hat;
    KMS_output.cU_EAM = c;
    KMS_output.CVU_EAM = CV;
    KMS_output.EIU_EAM = EI;
    KMS_output.flagU_EAM = flagU_opt;
    KMS_output.timeU_EAM = toc(tU_EAM);
end
if ((strcmp('one-sided-LB',type) == 1 || strcmp('two-sided',type) ==1)) && flag_feas == 1
    % Find minimum
    tL_EAM = tic;
     if KMSoptions.FeasAll
        [thetaL_hat,thetaL_optbound,c,CV,EI,flagL_opt]  = KMS_3_EAM(-p,-1,theta_feas,theta_init,c_init,CV_init,maxviol_init,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
     else
        [thetaL_hat,thetaL_optbound,c,CV,EI,flagL_opt]  = KMS_3_EAM(-p,-1,theta_feas(end,:),theta_init,c_init,CV_init,maxviol_init,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
     end
 
    % Save additional output to the KMS_output structure
    KMS_output.thetaL_EAM = thetaL_hat;
    KMS_output.cL_EAM = c;
    KMS_output.CVL_EAM = CV;
    KMS_output.EIL_EAM = EI;
    KMS_output.flagL_EAM = flagL_opt;
    KMS_output.timeL_EAM = toc(tL_EAM);
end

%% Directly Solve for Confidence Intervals
% In order to prove that the EAM algorithm recovered the true upper/lower
% bounds, we can solve for the upper/lower bounds directly using fmincon.  
% This method is not recommended in general as it is very slow.
if KMSoptions.direct_solve == 1
    options_fmincon.SpecifyConstraintGradient = false;

    if (strcmp('one-sided-UB',type) == 1 || strcmp('two-sided',type) ==1)
        % Find maximum
        tU_true = tic;        
        obj_true =  @(theta)KMS_41_TrueObjective(theta,p);
        constraint_true =  @(theta)KMS_42_TrueConstraint(theta,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
        problem = createOptimProblem('fmincon','x0',theta_0,...
        'objective',obj_true,'lb',LB_theta,'ub',UB_theta,...
        'nonlcon',constraint_true,'options',options_fmincon,'Aineq', A_theta, 'bineq', b_theta);
        [thetaU_hat_true,thetaU_optbound_true,flagU_opt_true,~,~] = run(options_multistart,problem,multistart_num);
        
        % Save output to KMS_output
        KMS_output.thetaU_optbound_true = -thetaU_optbound_true;  
        KMS_output.thetaU_true = thetaU_hat_true;
        KMS_output.timeU_true = toc(tU_true);
        
    end
    if (strcmp('one-sided-LB',type) == 1 || strcmp('two-sided',type) ==1)
        % Find minimum
        tL_true = tic;
        obj_true =  @(theta)KMS_41_TrueObjective(theta,-p);
        constraint_true =  @(theta)KMS_42_TrueConstraint(theta,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
        
        problem = createOptimProblem('fmincon','x0',theta_0,...
        'objective',obj_true,'lb',LB_theta,'ub',UB_theta,...
        'nonlcon',constraint_true,'options',options_fmincon,'Aineq', A_theta, 'bineq', b_theta);
        [thetaL_hat_true,thetaL_optbound_true,flagL_opt_true,~,~] = run(options_multistart,problem,multistart_num);
        
        % Save output to KMS_output
        KMS_output.thetaL_optbound_true = thetaL_optbound_true;
        KMS_output.thetaL_true = thetaL_hat_true;
        KMS_output.timeL_true = toc(tL_true);
    end
end

%% Output confidence interval
% We output confidence interval [thetaL,thetaU].  Note that if one of the
% bounds failed to converge, we output a warning message.  
%
% NOTE: If the one-sided LB/UB is requested, then we output the assumed
% bound on the parameter space
if  strcmp('two-sided',type) ==1
    thetaU = thetaU_optbound;
    thetaL = thetaL_optbound;  
    if flagU_opt ~= 1
        warning('Failed to find upper bound on the confidence interval -- confidence interval may be biased')
    end
    if flagL_opt ~= 1
        warning('Failed to find lower bound on the confidence interval -- confidence interval may be biased')
    end
elseif strcmp('one-sided-UB',type) == 1
    thetaU = thetaU_optbound;
    thetaL = CI_lo;
    if flagU_opt ~= 1
        warning('Failed to find upper bound on the confidence interval -- confidence interval may be biased')
    end
elseif strcmp('one-sided-LB',type) == 1
    thetaU = CI_hi;
    thetaL = thetaL_optbound;
    if flagL_opt ~= 1
        warning('Failed to find upper bound on the confidence interval -- confidence interval may be biasede')
    end
end
KMS_confidence_interval = [thetaL,thetaU];
end
