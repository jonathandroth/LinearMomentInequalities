function [g_ineq,g_eq] = moments_theta(theta,J1,J2,KMSoptions)
%% USER-SPECIFIED FUNCTION: Moment function that depends on parameter only
% The moment functions are in the form
%
%       E_P[m(W_i,theta)] = E_P[f(W_i)] + g(theta)
%
% where
%
%       E_P[m_j(W_i,theta)] <= 0 for j = 1,...,J1
%       E_P[m_j(W_i,theta)] = 0  for j = J1+1,...,J1+J2
%
% This function computes the model-implied component moment (in)equalities,
% g_j(theta).
%
% The user needs to specify this function.  An example is given below.
%
% INPUT:
%   theta         d-by-1 vector of parameters
%
%   J1            Integer number of moment inequalities
%
%   J2.           Integer number of moment equalities
%
%   KMSoptions.   This is a structure of additional inputs.  The user can
%                 add parameters to KMSoptions, say KMSoptions.params,
%                 and call KMSoptions.params in the user-specified
%                 functions.
%                 For example, in the 2x2 entry game, we include the
%                 support for the covariates and the probability that
%                 a particular support point occurs.
%
% OUTPUT:
%   g_ineq    J1-by-1 vector of moment inequalities g_j(theta) for j=1,...,J1
%
%   g_eq      2*J2-by-1 vector of moment inequalities g_j(theta) for j=1,...,J2
%             Note that we re-write the moment equalities as two moment
%             inequalities.  Thus, we have
%             f_eq = [g(theta) ; - g(theta)], where g(theta) is the
%             vector of moment equalities.
%
% NOTE: User does not need to specify J1 nor J2 here, since it is already
% computed in moments_w(W,KMSoptions).
%
% Below is a list of examples of moment functions.

% Select DGP
DGP = KMSoptions.DGP;

if DGP == -1
   %If DGP == -1, we assume that KMS options has a parameter X that is
   %saved, and the moments are of the form E[Y - X theta] >= 0
   %This function returns g(theta) = -X* theta
   g_ineq = -KMSoptions.X * theta;
   g_eq = [];
end


if DGP == 0
   g_ineq = [theta(1) + theta(2); theta(1) - theta(2); theta(1)];
   g_eq = [];
end

if DGP == 1
    %% Example:  Rotated cube (DGP-1)
    % Extract parameters
    theta_1 = theta(1);
    theta_2 = theta(2);
    
    % Define moment inequalities
    g_ineq = [theta_1 + theta_2 ; -theta_1 + theta_2 ; theta_1 - theta_2 ; -theta_1 - theta_2];
    
    % There are no moment equalities, so output the empty set.
    g_eq = [];

elseif DGP == 2
    %% Example:  Rotated cube (DGP-2)
    % Extract datasize and parameters
    n = KMSoptions.data_size;
    theta_1 = theta(1);
    theta_2 = theta(2);
      
    % Define moment inequalities
    g_ineq = [theta_1*sqrt(n) + theta_2 ; -theta_1*sqrt(n) + theta_2 ; theta_1*sqrt(n) - theta_2 ; -theta_1*sqrt(n) - theta_2];
    
    % There are no moment equalities, so output the empty set.
    g_eq = [];

elseif DGP == 2
    %% Example:  Rotated cube (DGP-3)
    % Extract datasize and parameters
    theta_1 = theta(1);
    theta_2 = theta(2);
    
    % Define moment inequalities
    g_ineq = [theta_1 + theta_2 ; -theta_1 + theta_2 ; theta_1 - theta_2 ; -theta_1 - theta_2];
    
    % There are no moment equalities, so output the empty set.
    g_eq = [];
    
elseif DGP == 4
    %% Example:  Rotated cube (DGP-1)
      % Extract datasize and parameters
    theta_1 = theta(1);
    theta_2 = theta(2);
    
    % Define moment inequalities
    g_ineq = [theta_1 + theta_2 ; -theta_1 + theta_2 ; theta_1 - theta_2 ; -theta_1 - theta_2;
              theta_1 + theta_2 ; -theta_1 + theta_2 ; theta_1 - theta_2 ; -theta_1 - theta_2];
    
    % There are no moment equalities, so output the empty set.
    g_eq = [];

elseif DGP == 5 || DGP == 6
    %% Example: 2-by-2 Entry Game
    % (See Pg 15)
    % Parameter vector is
    % theta =
    % (beta^1_0,            % Constant coefficient for firm 1
    % beta^1_1,             % Linear coefficient for firm 1
    % beta^2_0,             % Constant coefficient for firm 2
    % beta^2_1,             % Linear coefficient for firm 2
    % delta^1_0,            % Constant competition effect for firm 1
    % delta^1_1,            % Linear competition effect for firm 1
    % delta^2_0,            % Constant competition effect for firm 2
    % delta^2_1)            % Linear competition effect for firm 2
    
    % PRELIMINARY
    suppX = KMSoptions.suppX;                 % Support for covariates [x1,x2], which is [-1 -1; -1 1 ;  1 -1 ; 1  1 ];
    psuppX = KMSoptions.psuppX;               % Probability of support point occuring [0.1;0.2;0.3;0.4];
    dim_suppX = size(suppX,1);                % Number of support points
    g_ineq = zeros(J1,1);                     % Preset moment inequalities and
    g_eq = zeros(J2,1);                       % moment equalities

    % Extract parameters (easier to read)
    beta1  = theta(1:2);
    beta2  = theta(3:4);
    delta1 = theta(5:6);
    delta2 = theta(7:8);

    % MOMENT COMPUTATION
    % For each point in the support, we get the theoretical frequencies of
    % (Y1=y1,Y2=y2,X1=x1,X2=x2).
    % Potential outcomes are:
    % (y1,y2) = (0,0)  (both firms do not enter)
    % (y1,y2) = (1,1)  (both firms enter)
    % (y1,y2) = (1,0)  (only firm 1 enters)
    % (y1,y2) = (0,1)  (Only firm 2 enters)
    for ii = 1:dim_suppX
        % Pick support point (include constant)
        x1 = [1,suppX(ii,2)];
        x2 = [1,suppX(ii,4)];

        % Probability of this support point occuring
        pX= psuppX(ii);

        % Moment inequalities for entry game (See Pg 34, eq 5.3)
        g_ineq((ii-1)*2 + 1,1) = normcdf(-x1*(beta1+delta1))*(1-normcdf(-x2*beta2))*pX;

        % Moment inequalities for entry game (See Pg 34, eq 5.4)
        g_ineq((ii-1)*2 + 2,1) = ...
            -normcdf(-x1*(beta1+delta1))*(1-normcdf(-x2*(beta2+delta2)))*pX...
            -normcdf(-x1*beta1)*(normcdf(-x2*(beta2+delta2))-normcdf(-x2*beta2))*pX;

        % Moment equalities for entry game (See Pg 34, eq 5.1)
        g_eq((ii-1)*2 + 1,1)   = normcdf(-x1*beta1)*normcdf(-x2*beta2)*pX;

        % Moment equalities for entry game (See Pg 34, eq 5.2)
        g_eq((ii-1)*2 + 2,1) =(1-normcdf(-x1*(beta1+delta1)))*(1-normcdf(-x2*(beta2+delta2)))*pX;

    end
    % Concatenate the positive and negative of g_eq to transform into moment
    % inequalities
    g_eq = [g_eq ;-g_eq];

    % The way that the moment functions are written, m = f + g, g enters
    % additively.  The moment functions in this example have the incorrect
    % sign.
    g_ineq  = -g_ineq;
    g_eq    = -g_eq;

    % NOTE: In general the user needs to be careful in how s/he inputs the
    % moment conditions.  The program assumes that m = f+g.
    
elseif DGP == 7
    %% Example: 2-by-2 Entry Game With Correlated Errors
    % (See Pg 15)
    % Parameter vector is
    % theta =
    % (beta^1_0,            % Constant coefficient for firm 1
    % beta^1_1,             % Linear coefficient for firm 1
    % beta^2_0,             % Constant coefficient for firm 2
    % beta^2_1,             % Linear coefficient for firm 2
    % delta^1_0,            % Constant competition effect for firm 1
    % delta^1_1,            % Linear competition effect for firm 1
    % delta^2_0,            % Constant competition effect for firm 2
    % delta^2_1)            % Linear competition effect for firm 2
    % rho                   % Correlation in Bivariate Normal Distribution

    % PRELIMINARY
    suppX = KMSoptions.suppX;                   % Support for covariates [x1,x2], which is [-1 -1; -1 1 ;  1 -1 ; 1  1 ];
    psuppX = KMSoptions.psuppX;                 % Probability of support point occuring [0.1;0.2;0.3;0.4];
    dim_suppX = size(suppX,1);                  % Number of support points
    g_ineq = zeros(J1,1);                     % Preset moment inequalities and
    g_eq = zeros(J2,1);                       % moment equalities

    % Extract parameters (easier to read)
    beta1  = theta(1:2);
    beta2  = theta(3:4);
    delta1 = theta(5:6);
    delta2 = theta(7:8);
    rho    = theta(9);
    
    % MOMENT COMPUTATION
    % For each point in the support, we get the theoretical frequencies of
    % (Y1=y1,Y2=y2,X1=x1,X2=x2).
    % Potential outcomes are:
    % (y1,y2) = (0,0)  (both firms do not enter)
    % (y1,y2) = (1,1)  (both firms enter)
    % (y1,y2) = (1,0)  (only firm 1 enters)
    % (y1,y2) = (0,1)  (Only firm 2 enters)
    for ii = 1:dim_suppX
        % Pick support point (include constant)
        x1 = [1,suppX(ii,2)];
        x2 = [1,suppX(ii,4)];

        % Probability of this support point occuring
        pX= psuppX(ii);

        % Moment inequalities for entry game (See Pg 34, eq 5.3)
        g_ineq((ii-1)*2 + 1,1) = mexBVNcdf([-x1*(beta1+delta1),x2*beta2],[0;0],[1,-rho;-rho,1])*pX;
        
        % Moment inequalities for entry game (See Pg 34, eq 5.4)
        g_ineq((ii-1)*2 + 2,1) = ...
        -mexBVNcdf([-x1*(beta1+delta1),x2*(beta2+delta2)],[0;0],[1,-rho;-rho,1])*pX...
        -mexBVNcdf([-x1*beta1,-x2*(beta2+delta2)],[0;0],[1,rho;rho,1])*pX ...
        +mexBVNcdf([-x1*beta1,-x2*beta2],[0;0],[1,rho;rho,1])*pX;
    
        % Moment equalities for entry game (See Pg 34, eq 5.1)
        g_eq((ii-1)*2 + 1,1) = mexBVNcdf([-x1*beta1,-x2*beta2],[0;0],[1,rho;rho,1])*pX;
    
        % Moment equalities for entry game (See Pg 34, eq 5.2)
        g_eq((ii-1)*2 + 2,1) = mexBVNcdf([x1*(beta1+delta1),x2*(beta2+delta2)],[0;0],[1,rho;rho,1])*pX;
    end
    % Concatenate the positive and negative of g_eq to transform into moment
    % inequalities
    g_eq = [g_eq ;-g_eq];

    % The way that the moment functions are written, m = f + g, g enters
    % additively.  The moment functions in this example have the incorrect
    % sign.
    g_ineq  = -g_ineq;
    g_eq    = -g_eq;

    % NOTE: In general the user needs to be careful in how s/he inputs the
    % moment conditions.  The program assumes that m = f+g.
    
elseif DGP == 8   
    %% BCS Simulation
    % PRELIMINARY
    psuppX = KMSoptions.psuppX;                                             % Prob. of X_q=1 for q=1,2,3.
    dX = KMSoptions.dX;                                                     % Number of covariates
    g_ineq = zeros(J1,1);                                                   % Preset moment inequalities and 
    g_eq = zeros(J2,1);                                                     % moment equalities

    % Extract parameters (easier to read)
    theta1 = theta(1);
    theta2 = theta(2);
    beta1 = theta(3);
    beta2 = theta(4);
    beta3 = theta(5);
    beta = [0; beta1; beta2; beta3];

    % MOMENT COMPUTATION
    % For each point in the support, we get the theoretical frequencies (or bounds) of
    % (Y1=y1,Y2=y2,Xq=1).

    for ii = 1:dX
        g_ineq((ii-1)*2 + 1,1) = (theta2-beta(ii))*(1-theta1+beta(ii))*psuppX(ii);
        g_ineq((ii-1)*2 + 2,1) = -(theta2-beta(ii))*psuppX(ii);
        g_eq(ii,1)           = -(1-theta1+beta(ii))*(1-theta2+beta(ii))*psuppX(ii);
    end

    % Concatenate the positive and negative of g_eq to transform into moment
    % inequalities
    g_eq = [g_eq ;-g_eq];

    % NOTE: In general the user needs to be careful in how s/he inputs the
    % moment conditions.  The program assumes that m = f+g.  
end
end
