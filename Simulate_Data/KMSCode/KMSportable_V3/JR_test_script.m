%% Working Example: DGP6
% This script generates one data set for DGP6 with n=4000 observations.

method      = 'AS';    % Method - either AS or KMS
DGP         = 0;        % DGP = 0 (My added example)                     
alpha       = 0.05;     % Significance level
component   = 1;        % Component of theta to build confidence interval around
n           = 30;     % Sample size
KMSoptions  = KMSoptions();

KMSoptions.CVXGEN = 0;
KMSoptions.CVXGEN_name = [];


%% Extract/Save Information to KMSoptions, and set seed
KMSoptions.DGP = DGP;
KMSoptions.n = n;
KMSoptions.component = component;
seed = KMSoptions.seed;
B    = KMSoptions.B;            % Bootstraps
stream = RandStream('mlfg6331_64','Seed',seed);
RandStream.setGlobalStream(stream)

%% Parameters
type = 'two-sided';       % Two-sided or one sided test?  Set to 'one-sided-UB' or 'one-sided-LB' or 'two-sided'
kappa = NaN;              % Default kappa function
phi   = NaN;              % Default GMS function      

    
%% Run KMS

k = 3;
n= 500;
W = normrnd(0,1,n,k);

[f_ineq,f_eq,f_ineq_keep,f_eq_keep,paired_mom,J1,J2,J3] = moments_w(W,KMSoptions);
[g_ineq,g_eq] = moments_theta([1;-2],J1,J2,KMSoptions)
[f_stdev_ineq,f_stdev_eq] = moments_stdev(W,f_ineq,f_eq,J1,J2,KMSoptions)
[Dg_ineq,Dg_eq] = moments_gradient([15;-2],J1,J2,KMSoptions)

LB_theta = [-10;-10];
UB_theta = [10;10];
alpha = 0.05;
CVXGEN_name = [];
theta_0 = [0;0];
component = 1;
p                   = zeros(size(theta_0,1),1);                            % Projection direction
p(component)        = 1;

%KMSoptions.parallel = 0; %turn off parallelization for debugging.
KMSoptions.S = 0; %the default S=[] seems to cause dimension issues? 
KMSoptions.B = 1000; %set bootstrap draws

tic()
[KMS_confidence_interval,KMS_output] = KMS_0_Main(W,theta_0,...
    p,[],LB_theta,UB_theta,[],[],alpha,type,method,kappa,phi,CVXGEN_name,KMSoptions);
toc()

%% Save an X in KMS options and run the KMS code with dgp == -1, corresponding with E[Y + X delta]>=0
KMSoptions.X = [-1,0;1,0;0,0];
KMSoptions.DGP = -1;
p = [1;0];
tic()
[KMS_confidence_interval,KMS_output] = KMS_0_Main(W,theta_0,...
    p,[],LB_theta,UB_theta,[],[],alpha,type,method,kappa,phi,CVXGEN_name,KMSoptions);
toc()


%% Try out the AS function
m = 5; %number of moments
k = 3; %number of parameters
Y = normrnd(0,1,n,m);
Y_T = mean(Y)'* sqrt(n);
X_T = normrnd(0,1,m,k);
projected_AS(Y_T, X_T, n, eye(m), [1;0;0])

y_T = [1;-1;0]*sqrt(n);
X_T = [ [-1;1;0] , [0;0;0] ];
projected_AS(y_T, X_T, n, eye(3), [1;0])