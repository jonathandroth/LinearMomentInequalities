%% Code description
% This code runs simulations using the KMS_0_Main code.
% This can be used as example code for either the games or rotated box.  
%
% Key Parameters to set below:
%   Component       The component of the paramter vector that we build
%                   projected confidence sets for.
%  
%   DGP             Equal to 1,...,7.  DGP  = 1,2,3,4 are rotated-box
%                   examples, DGP = 5,6 are games examples, 
%                   (point-identified case and partially-identified case),
%                   and DGP=7 is the BCS example. 
%   
%   alpha           Significance level
%
%   method          Equal to 'KMS' or 'AS', depending on the method to
%                   build the confidence set.  'AS' tends to overinflate
%                   the confidence set
%
%   n               Sample size
%
%   name            Name of file to save to results folder
%
%   KMSoptions      Structure of options set by user.
clear
clc

%% Key Parameters:
method      = 'KMS';    % Method - either AS or KMS
DGP         = 1;        % DGP = 1,2,3,4,5,6,7,8.  DGP = 5,6,7 are Games examples.  
                        % DGP=8 is the BCS example.  If DGP = 8 and BCS = 1,
                        % then also calls computation for BCS confidence interval 
KMS         = 8;        % Set equal to 1 to run KMS simulations.  
                        % Set DGP = 8 and KMS = 0 to run BCS. 
alpha       = 0.15;     % Significance level
component   = 1;        % Component of theta to build confidence interval around
n           = 4000;     % Sample size
Nmc         = 1;        % Number of Monte Carlos      
sim_lo      = 1;        % BCS code is very slow.  We split the job between many computers.
sim_hi      = Nmc;      % We run simulations mm = sim_lo ... sim_hi.
                        % Default is sim_lo = 1 and sim_hi= Nmc
name        = strcat(method,'_DGP=',num2str(DGP),'_coverage=',num2str(100*(1-alpha)),'_component=',num2str(component),'_samplesize',num2str(n),'_numMC',num2str(Nmc),'_simnums=',num2str(sim_lo),'to',num2str(sim_hi),'_');
KMSoptions  = KMSoptions();

%% Extract/Save Information to KMSoptions, and set seed
KMSoptions.DGP = DGP;
KMSoptions.n = n;
KMSoptions.component = component;
seed = KMSoptions.seed;
B    = KMSoptions.B;                                                        % Bootstraps
stream = RandStream('mlfg6331_64','Seed',seed);
RandStream.setGlobalStream(stream)

%% Parameters
type = 'two-sided';         % Two-sided or one sided test?  Set to 'one-sided-UB' or 'one-sided-LB' or 'two-sided'
kappa =NaN;                 % Default kappa function
phi   = NaN;                % Default GMS function      

%% Parameters that depend on DGP
if DGP == 1 || DGP == 2 || DGP == 3
    theta_0 = [0;-1];               % Initial parameter vector
    p = zeros(size(theta_0,1),1);   % Projection direction
    p(component) = 1;               % Projection direction
    KMSoptions.S =  0;              % Rho Polytope Constraints 
    LB_theta = [-2; -4];            % Lower bound on parameter space
    UB_theta = [2;2];               % Upper on parameter space
    A_theta = [];                   % No polytope constraints
    b_theta = [];                   % No polytope constraints
    CVXGEN_name = 'csolve_DGP1';    % CVXGEN file name
    
elseif DGP == 4
    theta_0 = [0;-1];               % Initial parameter vector
    p = zeros(size(theta_0,1),1);   % Projection direction
    p(component) = 1;               % Projection direction
    KMSoptions.S =  0;              % Rho Polytope Constraints
    LB_theta = [-2; -4];            % Lower bound on parameter space
    UB_theta = [2;2];               % Upper bound on parameter space
    A_theta = [];                   % No polytope constraints
    b_theta = [];                   % No polytope constraints
    CVXGEN_name = 'csolve_DGP4';    % CVXGEN file name
    
elseif DGP == 5 
    theta_true = [0.5 ; 0.25; 0.5 ; 0.25; -1; -1 ; -1 ;-1];                % True parameter vector
    LB_theta = [-2*ones(4,1); -4*ones(4,1)];                               % Lower bound on parameter space
    UB_theta = [2*ones(4,1); zeros(4,1)];                                  % Upper bound on parameter space
    theta_0 = 0.5*UB_theta + 0.5*LB_theta;                                 % Set initial parameter vector to something feasible.
    p = zeros(size(theta_0,1),1);                                          % Projection direction
    p(component) = 1;                                                      % Projection direction
    KMSoptions.S = 0;                                                      % Rho Polytope Constraints
    A_theta = [];                                                          % No polytope constraints
    b_theta = [];                                                          % No polytope constraints
    suppX = [1 -1 1 -1 ; 1 -1 1 1 ; 1 1 1 -1 ; 1 1 1 1];                   % Support for X
    psuppX = [0.1;0.2;0.3;0.4];                                            % Prob of support point occuring
    KMSoptions.suppX = suppX;
    KMSoptions.psuppX = psuppX;
    selp = 0.5;                                                            % Selection mechanism for non-unique equilibria Pr(Y = (0,1))
    KMSoptions.selp = selp;
    CVXGEN_name = 'csolve_DGP5';                                           % CVXGEN file name
    
elseif DGP == 6 
    theta_true = [0.5 ; 0.25; 0.5 ; 0.25; -1; -0.75 ; -1 ;-0.75];          % True parameter vector
    LB_theta = [-2*ones(4,1); -4*ones(4,1)];                               % Lower bound on parameter space
    UB_theta = [2*ones(4,1); zeros(4,1)];                                  % Upper bound on parameter space
    %theta_0 = theta_true;                                                 % Initial parameter vector (not feasible in practice)
    theta_0 = 0.5*UB_theta + 0.5*LB_theta;                                 % Set initial parameter vector to something feasible.
    p = zeros(size(theta_0,1),1);                                          % Projection direction
    p(component) = 1;                                                      % Projection direction
    KMSoptions.S =  0;                                                     % Rho Polytope Constraints 
    A_theta = [];
    b_theta = [];
    suppX = [1 -1 1 -1 ; 1 -1 1 1 ; 1 1 1 -1 ; 1 1 1 1];                   % Support for X
    psuppX = [0.1;0.2;0.3;0.4];                                            % Prob of support point occuring
    KMSoptions.suppX = suppX;
    KMSoptions.psuppX = psuppX;
    selp = 0.5;                                                             % Selection mechanism for non-unique equilibria Pr(Y = (1,0))
    KMSoptions.selp = selp;
    CVXGEN_name = 'csolve_DGP5';                                            % CVXGEN file name

elseif DGP == 7
    theta_true  = [0.5 ; 0.25; 0.5 ; 0.25; -1; -0.75 ; -1 ;-0.75 ; 0.5];   % True parameter vector
    LB_theta = [-2*ones(4,1); -4*ones(4,1) ; 0];                           % Lower bound on parameter space
    UB_theta = [2*ones(4,1); zeros(4,1); 0.95];                            % Upper bound on parameter space
    theta_0 = 0.5*UB_theta + 0.5*LB_theta;                                 % Set initial parameter vector to something feasible.
    p = zeros(size(theta_0,1),1);                                          % Projection direction
    p(component) = 1;
    KMSoptions.S =  0;                                                     % Rho Polytope Constraints
    A_theta = [];
    b_theta = [];
    suppX = [1 -1 1 -1 ; 1 -1 1 1 ; 1 1 1 -1 ; 1 1 1 1];                   % Support for X
    psuppX = [0.1;0.2;0.3;0.4];                                            % Prob of support point occuring
    KMSoptions.suppX = suppX;
    KMSoptions.psuppX = psuppX;
    selp = 0.5;                                                            % Selection mechanism for non-unique equilibria Pr(Y = (1,0))
    KMSoptions.selp = selp;
    CVXGEN_name = 'csolve_DGP7';    
    
elseif DGP == 8   
    theta_true  = [0.4 ; 0.6 ;0.1  ;0.2  ;0.3];                            % True parameter vector
    dim_p       = size(theta_true,1);
    p = zeros(size(theta_true,1),1);                                       % Projection direction
    p(component) = 1;
    KMSoptions.S =  13;                                                    % Rho Polytope Constraints
    % Param space is theta_1,theta_2 in [0,1], theta_k in
    % [0,min(theta_1,theta_2)] for k = 1,2,3.
    LB_theta    = zeros(dim_p,1);
    UB_theta    = ones(dim_p,1);
    A_theta     = [-1 0 1 0 0 ; 0 -1 1 0 0  ; -1 0 0 1 0 ; 0 -1 0 1 0  ;  -1 0 0 0 1 ; 0 -1 0 0 1 ]; 
    b_theta     = zeros(6,1);
    % We randomly select theta_0 from the parameter space
    theta_0     =0.5*LB_theta + 0.5*UB_theta;   
    dX = size(theta_true,1)-1;                                             % dimension of X.
    psuppX = ones(dX,1)/dX;                                                % P(X=x), X is discrete uniform.
    KMSoptions.dX = dX;
    KMSoptions.psuppX = psuppX;
    selp = 0.6;                                                            % prob of P(A_1=0,A_2=1) when there is multiplicity.
    KMSoptions.selp = selp;                                                % NOTE: THIS IS THE OPPOSITE of DGP6,DGP5.
    CVXGEN_name = 'csolve_DGP8';                                           % CVXGEN file name
end

%% Generate data
if DGP == 1 || DGP == 2 || DGP == 3 
    % DATA FOR ROTATED BOX:
    % u1 and u1 are random normal preference shocks that determine payoff in
    % states where the firm decides to enter
    Z = randn(n,4,Nmc);
elseif DGP == 4
    % DATA FOR ROTATED BOX:
    Z = randn(n,8,Nmc);
    Z(:,[5,7],:)=2.*Z(:,[5,7],:);
    Z(:,[6,8],:)=3.*Z(:,[6,8],:);
elseif DGP == 5 || DGP == 6
    % DATA FOR GAME:
    % u1 and u1 are random normal preference shocks that determine payoff in
    % states where the firm decides to enter
    u1 = randn(n,Nmc);
    u2 = randn(n,Nmc);

    % v is a latent binary variable that determines which equilibrium are
    % selected in the multi-equilibrium case
    v  = (rand(n,Nmc)>selp);

    % X = (x1, x2) are characteristic of firms, which are distributed on suppX.  It
    % is assumed that x1, x2 are binary, hence we use the gendist function to
    % generate random variables from a discrete probability distribution.
    X = zeros(n,4,Nmc);
    temp = KMS_AUX3_gendist(psuppX.',n,Nmc);
    for mm = 1:Nmc
        X(:,:,mm) = suppX(temp(:,mm),:);
    end

    % Determine outcomes Y = (y1,y2)
    % Potential outcomes are:
    % (y1,y2) = (0,0)  (both firms do not enter)
    % (y1,y2) = (1,1)  (both firms enter)
    % (y1,y2) = (1,0)  (only firm 1 enters)
    % (y1,y2) = (0,1)  (Only firm 2 enters)
    Y = zeros(n,2,Nmc);

    % Extract parameters (easier to read)
    beta1 = theta_true(1:2);
    beta2 = theta_true(3:4);
    delta1 = theta_true(5:6);
    delta2 = theta_true(7:8);
    
    for mm = 1:Nmc
        for ii = 1:n
            % For each market and each simulation run, we determine the outcome
            x1 = X(ii,1:2,mm);                           % Covariate
            x2 = X(ii,3:4,mm);                           % Covariate

            if u1(ii,mm) <= -x1*beta1 && u2(ii,mm) <= -x2*beta2
                % (1) unique equlibrium is (0,0)
                Y(ii,:,mm) = [0,0];

            elseif u1(ii,mm) >= -x1*(beta1+delta1) && u2(ii,mm) >= -x2*(beta2+delta2)
                % (2) unique equlibirum is (1,1)
                Y(ii,:,mm) = [1,1];

            elseif (u1(ii,mm) <= -x1*(beta1+delta1) && u2(ii,mm) >= -x2*(beta2+delta2)) || ...
                    (u1(ii,mm) <= -x1*beta1 && u2(ii,mm) <= -x2*(beta2+delta2) && u2(ii,mm)>= -x2*beta2)
                % (3) unique equlibirum is (0,1)
                Y(ii,:,mm) = [0,1];

            elseif u1(ii,mm) >= -x1*beta1 && u1(ii,mm) < -x1*(beta1 + delta1) ...
                    && u2(ii,mm) >= -x2*beta2 && u2(ii,mm) <= -x2*(beta2+delta2)
                % (4) multiple equilibria region, select (0,1) with probability
                % selp
                if v(ii,mm) == 1
                    Y(ii,:,mm) = [1,0];
                else
                    Y(ii,:,mm) = [0,1];
                end

            else
                % (5) unique equilibrium is (1,0)
                Y(ii,:,mm) = [1,0];
            end
        end
    end
    Z = [Y, X];
    
elseif DGP == 7
    % DATA FOR GAME:
    % u1 and u1 are random normal preference shocks that determine payoff in
    % states where the firm decides to enter.  In this DGP u1 and u2 are
    % normal random variables with correlation coefficient equal to rho
    % (last component of theta)
    mu = [0,0];
    sigma = [1,theta_true(end);theta_true(end),1];
    U=mvnrnd(mu,sigma,n*Nmc);
    u1=reshape(U(:,1),[n,Nmc]);
    u2=reshape(U(:,2),[n,Nmc]);

    % v is a latent binary variable that determines which equilibrium are
    % selected in the multi-equilibrium case
    v  = (rand(n,Nmc)>selp);

    % X = (x1, x2) are characteristic of firms, which are distributed on suppX.  It
    % is assumed that x1, x2 are binary, hence we use the gendist function to
    % generate random variables from a discrete probability distribution.
    X = zeros(n,4,Nmc);
    temp = KMS_AUX3_gendist(psuppX.',n,Nmc);
    for mm = 1:Nmc
        X(:,:,mm) = suppX(temp(:,mm),:);
    end

    % Determine outcomes Y = (y1,y2)
    % Potential outcomes are:
    % (y1,y2) = (0,0)  (both firms do not enter)
    % (y1,y2) = (1,1)  (both firms enter)
    % (y1,y2) = (1,0)  (only firm 1 enters)
    % (y1,y2) = (0,1)  (Only firm 2 enters)
    Y = zeros(n,2,Nmc);

    % Extract parameters (easier to read)
    beta1 = theta_true(1:2);
    beta2 = theta_true(3:4);
    delta1 = theta_true(5:6);
    delta2 = theta_true(7:8);
    
    for mm = 1:Nmc
        for ii = 1:n
            % For each market and each simulation run, we determine the outcome
            x1 = X(ii,1:2,mm);                           % Covariate
            x2 = X(ii,3:4,mm);                           % Covariate

            if u1(ii,mm) <= -x1*beta1 && u2(ii,mm) <= -x2*beta2
                % (1) unique equlibrium is (0,0)
                Y(ii,:,mm) = [0,0];

            elseif u1(ii,mm) >= -x1*(beta1+delta1) && u2(ii,mm) >= -x2*(beta2+delta2)
                % (2) unique equlibirum is (1,1)
                Y(ii,:,mm) = [1,1];

            elseif (u1(ii,mm) <= -x1*(beta1+delta1) && u2(ii,mm) >= -x2*(beta2+delta2)) || ...
                    (u1(ii,mm) <= -x1*beta1 && u2(ii,mm) <= -x2*(beta2+delta2) && u2(ii,mm)>= -x2*beta2)
                % (3) unique equlibirum is (0,1)
                Y(ii,:,mm) = [0,1];

            elseif u1(ii,mm) >= -x1*beta1 && u1(ii,mm) < -x1*(beta1 + delta1) ...
                    && u2(ii,mm) >= -x2*beta2 && u2(ii,mm) <= -x2*(beta2+delta2)
                % (4) multiple equilibria region, select (0,1) with probability
                % selp
                if v(ii,mm) == 1
                    Y(ii,:,mm) = [1,0];
                else
                    Y(ii,:,mm) = [0,1];
                end

            else
                % (5) unique equilibrium is (1,0)
                Y(ii,:,mm) = [1,0];
            end
        end
    end
    Z = [Y, X];
    
elseif DGP == 8   
    % DATA FOR BCS SIMULATION
    % Draw random variable
    data_BCS = zeros(n,2*dX,Nmc);
    Z = zeros(n,2*dX,Nmc);
    baseDatas = rand(n,4,Nmc); 
    for mm = 1:Nmc
        epsilons = baseDatas(:,1:2,mm);    % epsilon in the model in section 5;
        multiple = baseDatas(:,3,mm);      % determines how multiplicity is resolved;

        % X denotes the market type indicator;
        X = zeros(n,1);
        for j=1:dX
            X = X + (j-1)*(baseDatas(:,4,mm)>=(j-1)/dX).*(baseDatas(:,4,mm)<j/dX) ;
        end
        betas_true = [0;theta_true(3:end)]; % vector indicates beta_q for q=1,...,d_X

        % Initialize matrices that will contain both entry decision and market type
        dataP11 = zeros(n,dX);
        dataP10 = zeros(n,dX);

        dataP11_KMS = zeros(n,dX);
        dataP10_KMS = zeros(n,dX);

        for j=1:dX
            % Entry decision that indices {A_1=1,A_2=1}
            Entry11_aux = (epsilons(:,1) > theta_true(1)-betas_true(j)).*(epsilons(:,2) >  theta_true(2)-betas_true(j)) ;

            % Entry decision that indices {A_1=1,A_2=0}
            Entry10_aux = (epsilons(:,1) > theta_true(1)-betas_true(j)).*(epsilons(:,2) <= theta_true(2)-betas_true(j)) + (epsilons(:,1) <= theta_true(1)-betas_true(j)).*(epsilons(:,2) <= theta_true(2)-betas_true(j)).*(multiple>selp);

            % Data for BCS
            dataP11(:,j) = Entry11_aux.*(X == j-1)/psuppX(j);               % indicates {A_1=1,A_2=1} and {X=j-1}
            dataP10(:,j) = Entry10_aux.*(X == j-1)/psuppX(j);               % indicates {A_1=1,A_2=0} and {X=j-1}

            % KMS can work with unconditional moments
            dataP11_KMS(:,j) = Entry11_aux.*(X == j-1);  
            dataP10_KMS(:,j) = Entry10_aux.*(X == j-1);
        end
        data_BCS(:,:,mm)    = [dataP11, dataP10];
        Z(:,:,mm)           = [dataP11_KMS, dataP10_KMS];                  % Save Z to pass to KMS algorithm
    end
end

%% Compute population identification region
% This step is only relevant for simulation.  
if DGP == 5 || DGP == 6 || DGP == 7 || DGP == 8
    stream = RandStream('mlfg6331_64','Seed',seed);
    RandStream.setGlobalStream(stream)
    stream.Substream = B + B*10^3 + 2;
    addpath ./MVNorm
    Identification_region = KMS_5_identification_region(theta_true,theta_0,LB_theta,UB_theta,A_theta,b_theta,KMSoptions);
    KMSoptions.Identification_region = Identification_region;
    stream = RandStream('mlfg6331_64','Seed',seed);
    RandStream.setGlobalStream(stream)
    stream.Substream = B + B*10^3 + 3;
end

if KMS
    %% Run KMS
    warning('off')
    t1 = tic;
    for mm = sim_lo:sim_hi
        fprintf(['Starting Monte Carlo Number for ', method,'  %d \n'],mm)
        t2 = tic;
        W = Z(:,:,mm);
        [KMS_confidence_interval{mm},KMS_output{mm}] = KMS_0_Main(W,theta_0,...
            p,[],LB_theta,UB_theta,A_theta,b_theta,alpha,type,method,kappa,phi,CVXGEN_name,KMSoptions);
        KMS_output{mm}.totaltime = toc(t2);
    end
    totaltime_KMS = toc(t1);

    %% Cell to vector
    KMS_CI = nan(Nmc,2);
    for mm =  sim_lo:sim_hi
        KMS_CI(mm,:) = KMS_confidence_interval{mm};
    end

    %% Save  KMS
    date = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
    filename = strcat('Results/KMSresults_',name,date,'.mat');
    save(filename)
end

% Run BCS (if required)
if DGP == 8 && (component == 1 || component == 2) && ~KMS
    addpath ./BCS
    t1 = tic;
    for mm = sim_lo:sim_hi
        fprintf('Starting Monte Carlo Number for BCS %d \n',mm)
        t2 = tic;
        W = data_BCS(:,:,mm);
        [BCS_confidence_interval{mm},BCS_output{mm}] = BCS_Main(W,theta_true,component,selp,n,B,seed,alpha);
        BCS_output{mm}.totaltime = toc(t2);
    end
    totaltime_BCS = toc(t1);
    
    %% Save  BCS
    date = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
    filename = strcat('Results/BCSresults_',name,date,'.mat');
    save(filename)
end



