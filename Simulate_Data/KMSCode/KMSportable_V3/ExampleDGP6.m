%% Working Example: DGP6
% This script generates one data set for DGP6 with n=4000 observations.

method      = 'KMS';    % Method - either AS or KMS
DGP         = 6;        % DGP = 6 (Partially-Identified Games Example)                     
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

%% Parameters that depend on DGP
theta_true          = [0.5 ; 0.25; 0.5 ; 0.25; -1; -0.75 ; -1 ;-0.75];     % True parameter vector
LB_theta            = [-2*ones(4,1); -4*ones(4,1)];                        % Lower bound on parameter space
UB_theta            = [2*ones(4,1); zeros(4,1)];                           % Upper bound on parameter space
theta_0             = 0.5*UB_theta + 0.5*LB_theta;                         % Set initial parameter vector to something feasible.
p                   = zeros(size(theta_0,1),1);                            % Projection direction
p(component)        = 1;
KMSoptions.S        =  0;                                                  % Rho Polytope Constraints 
suppX               = [1 -1 1 -1 ; 1 -1 1 1 ; 1 1 1 -1 ; 1 1 1 1];         % Support for X
psuppX              = [0.1;0.2;0.3;0.4];                                   % Prob of support point occuring
KMSoptions.suppX    = suppX;
KMSoptions.psuppX   = psuppX;
selp                = 0.5;                                                 % Selection mechanism for non-unique equilibria Pr(Y = (1,0))
KMSoptions.selp     = selp;
%CVXGEN_name         = 'ExampleDGP6';                                       % CVXGEN file name
CVXGEN_name         = [];
%% Generate data
% DATA FOR GAME:
% u1 and u1 are random normal preference shocks that determine payoff in
% states where the firm decides to enter.  In this DGP u1 and u2 are
% normal random variables with correlation coefficient equal to rho
% (last component of theta)
mu = [0,0];
sigma = [1,theta_true(end);theta_true(end),1];
U=mvnrnd(mu,sigma,n);
u1=U(:,1);
u2=U(:,2);

% v is a latent binary variable that determines which equilibrium are
% selected in the multi-equilibrium case
v  = (rand(n,1)>selp);

% X = (x1, x2) are characteristic of firms, which are distributed on suppX.  It
% is assumed that x1, x2 are binary, hence we use the gendist function to
% generate random variables from a discrete probability distribution.
X = zeros(n,4);
temp = KMS_AUX3_gendist(psuppX.',n,1);
X = suppX(temp(:,1),:);

% Determine outcomes Y = (y1,y2)
% Potential outcomes are:
% (y1,y2) = (0,0)  (both firms do not enter)
% (y1,y2) = (1,1)  (both firms enter)
% (y1,y2) = (1,0)  (only firm 1 enters)
% (y1,y2) = (0,1)  (Only firm 2 enters)
Y = zeros(n,2);

% Extract parameters (easier to read)
beta1 = theta_true(1:2);
beta2 = theta_true(3:4);
delta1 = theta_true(5:6);
delta2 = theta_true(7:8);

for ii = 1:n
    % For each market and each simulation run, we determine the outcome
    x1 = X(ii,1:2);                           % Covariate
    x2 = X(ii,3:4);                           % Covariate

    if u1(ii) <= -x1*beta1 && u2(ii) <= -x2*beta2
        % (1) unique equlibrium is (0,0)
        Y(ii,:) = [0,0];

    elseif u1(ii) >= -x1*(beta1+delta1) && u2(ii) >= -x2*(beta2+delta2)
        % (2) unique equlibirum is (1,1)
        Y(ii,:) = [1,1];

    elseif (u1(ii) <= -x1*(beta1+delta1) && u2(ii) >= -x2*(beta2+delta2)) || ...
            (u1(ii) <= -x1*beta1 && u2(ii) <= -x2*(beta2+delta2) && u2(ii)>= -x2*beta2)
        % (3) unique equlibirum is (0,1)
        Y(ii,:) = [0,1];

    elseif u1(ii) >= -x1*beta1 && u1(ii) < -x1*(beta1 + delta1) ...
            && u2(ii) >= -x2*beta2 && u2(ii) <= -x2*(beta2+delta2)
        % (4) multiple equilibria region, select (0,1) with probability
        % selp
        if v(ii) == 1
            Y(ii,:) = [1,0];
        else
            Y(ii,:) = [0,1];
        end

    else
        % (5) unique equilibrium is (1,0)
        Y(ii,:) = [1,0];
    end
end
W = [Y, X];
    
%% Run KMS
[KMS_confidence_interval,KMS_output] = KMS_0_Main(W,theta_0,...
    p,[],LB_theta,UB_theta,[],[],alpha,type,method,kappa,phi,CVXGEN_name,KMSoptions);



