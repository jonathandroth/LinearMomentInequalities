function [Dg_ineq,Dg_eq] = moments_gradient(theta,J1,J2,KMSoptions)
%% USER-SPECIFIED FUNCTION: Gradient of the moment function
% The moment functions are in the form
%
%       E_P[m(W_i,theta)] = E_P[f(W_i)] + g(theta)
%
% where
%
%       E_P[m_j(W_i,theta)] <= 0 for j = 1,...,J1
%       E_P[m_j(W_i,theta)] = 0  for j = J1+1,...,J1+J2
%
% This function computes Dm(W_i,theta) = Dg(theta).
%
% NOTE: We used the separability assumption above, so that Dm does not
% depend on data.
%
% INPUT:
%   theta         dim_p-by-1 vector of parameters
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

% OUTPUT:
%   Dg_ineq     J1-by-dim_p vector of gradients of moment inequalities
%                 Dg_j(theta) for j = 1,...,J1
%
%   Dg_eq       2*J2-by-dim_p vector of gradients of moment equalities
%                 Dg_j(theta) for j = 1,...,J2
%
% Below is a list of examples of moment functions.

% Select DGP
DGP = KMSoptions.DGP;

if DGP == -1
   %If DGP == -1, we assume that KMS options has a parameter X that is
   %saved, and the moments are of the form E[Y - X theta] <= 0
   %This function returns the gradient of g(theta) = -X theta wrt theta, i.e. -X
   Dg_ineq = -KMSoptions.X;
   Dg_eq = [];
end


if DGP == 0
   Dg_ineq = [1,1; 1,-1; 1,0];
   Dg_eq = [];
end
if DGP == 1
    %% Example: Rotaed cube (DGP-1)
    Dg_ineq = [1 1 ; -1 1 ; 1 -1 ; -1 -1];
    Dg_eq   = [];

elseif DGP == 2
    %% Example: Rotaed cube (DGP-2)
    n = KMSoptions.data_size;
    Dg_ineq = [sqrt(n) 1 ; -sqrt(n) 1 ; sqrt(n) -1 ; -sqrt(n) -1];
    Dg_eq   = [];
    
elseif DGP == 3
    %% Example: Rotaed cube (DGP-3)
    Dg_ineq = [1 1 ; -1 1 ; 1 -1 ; -1 -1];
    Dg_eq   = [];

elseif DGP == 4
    %% Example: Rotaed cube (DGP-4)
    Dg_ineq = [1 1 ; -1 1 ; 1 -1 ; -1 -1;
               1 1 ; -1 1 ; 1 -1 ; -1 -1];
    Dg_eq   = [];

elseif DGP == 5 || DGP == 6 
    %% Example: 2-by-2 Entry Game
    % (See Pg 15)
    % Parameter vector is
    % theta =
    % (beta^1_0,            % Constant coefficient for firm 1
    % beta^1_1,             % Linear coefficient for firm 1
    % (beta^2_0,            % Constant coefficient for firm 2
    % beta^2_1,             % Linear coefficient for firm 2
    % delta^1_0,            % Constant competition effect for firm 1
    % delta^1_1,            % Linear competition effect for firm 1
    % delta^2_0,            % Constant competition effect for firm 2
    % delta^2_1)            % Linear competition effect for firm 2

    % PRELIMINARY
    dim_p = KMSoptions.dim_p;
    suppX = KMSoptions.suppX;                   % Support for covariates [x1,x2], which is [-1 -1; -1 1 ;  1 -1 ; 1  1 ];
    psuppX = KMSoptions.psuppX;                 % Probability of support point occuring [0.1;0.2;0.3;0.4];
    dim_suppX = size(suppX,1);                  % Number of support points
    Dg_ineq = zeros(J1,dim_p);                     % Preset moment inequalities and
    Dg_eq = zeros(J2,dim_p);                       % moment equalities

    % Extract parameters (easier to read)
    beta1 = theta(1:2);
    beta2 = theta(3:4);
    delta1 = theta(5:6);
    delta2 = theta(7:8);

    % GRADIENT COMPUTATION
    % For each point in the support, we get the theoretical gradient for each
    % potential point in the support
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

        % Probability of this support occuring
        pX= psuppX(ii);

        % Gradients for moments inequalities for entry game (See Pg 34, eq 5.3)
        Dg3b1 = x1.*normpdf(-x1*(beta1+delta1)) .*(1-normcdf(-x2*beta2))*pX;
        Dg3b2 = -x2.*normcdf(-x1*(beta1+delta1)).*normpdf(-x2*beta2)*pX;
        Dg3d1 = x1.*normpdf(-x1*(beta1+delta1)) .*(1-normcdf(-x2*beta2))*pX;
        Dg3d2 = zeros(1,2);

        % Gradients for moments inequalities for entry game (See Pg 34, eq 5.4)
        Dg4b1 = -x1.*normpdf(-x1*(beta1+delta1)).*(1-normcdf(-x2*(beta2+delta2)))*pX ...
               -x1.*normpdf(-x1*beta1).*(normcdf(-x2*(beta2+delta2))-normcdf(-x2*beta2))*pX;
        Dg4b2 = x2.*normcdf(-x1*(beta1+delta1)).*normpdf(-x2*(beta2+delta2))*pX ...
               +x2.*normcdf(-x1*beta1).*(-normpdf(-x2*(beta2+delta2))+normpdf(-x2*beta2))*pX;
        Dg4d1 = -x1.*normpdf(-x1*(beta1+delta1)).*(1-normcdf(-x2*(beta2+delta2)))*pX;
        Dg4d2 = x2.*normcdf(-x1*(beta1+delta1)).*normpdf(-x2*(beta2+delta2))*pX ...
              +x2.*normcdf(-x1*beta1).*(-normpdf(-x2*(beta2+delta2)))*pX;


        % Gradients for moment equalities for entry game (See Pg 34, eq 5.1)
        Dg1b1 = x1.*normpdf(-x1*beta1).*normcdf(-x2*beta2)*pX;
        Dg1b2 = x2.*normcdf(-x1*beta1).*normpdf(-x2*beta2)*pX;
        Dg1d1 = zeros(1,2);
        Dg1d2 = zeros(1,2);

        % Gradients for moment equalities for entry game (See Pg 34, eq 5.2)
        Dg2b1 = -x1.*normpdf(-x1*(beta1+delta1))     .*(1-normcdf(-x2*(beta2+delta2)))*pX;
        Dg2b2 = -x2.*(1-normcdf(-x1*(beta1+delta1))) .*normpdf(-x2*(beta2+delta2))*pX;
        Dg2d1 = -x1.*normpdf(-x1*(beta1+delta1))     .*(1-normcdf(-x2*(beta2+delta2)))*pX;
        Dg2d2 = -x2.*(1-normcdf(-x1*(beta1+delta1))) .*normpdf(-x2*(beta2+delta2))*pX;

        % Compile moments
        Dg_ineq((ii-1)*2 + 1,:) = [Dg3b1 Dg3b2 Dg3d1 Dg3d2];
        Dg_ineq((ii-1)*2 + 2,:) = [Dg4b1 Dg4b2 Dg4d1 Dg4d2];
        Dg_eq((ii-1)*2 + 1,:)   = [Dg1b1 Dg1b2 Dg1d1 Dg1d2];
        Dg_eq((ii-1)*2 + 2,:)   = [Dg2b1 Dg2b2 Dg2d1 Dg2d2];

        % NOTE: These gradients are indeed J1 or 2*J2-by-dim_p, since G3b1 ect
        % is 2-by-1 (x1 and x2 include constant).
    end
    % Concatenate the positive and negative of Dmom_eq to transform into moment
    % inequalities
    Dg_eq = [Dg_eq ; -Dg_eq];

    % The way that the moment functions are written, m = f + g, g enters
    % additively.  The moment functions in this example have the incorrect
    % sign.  However, the gradient functions have the correct sign!
    
elseif DGP == 7
     %% Example: 2-by-2 Entry Game With Correlated Errors
    % (See Pg 15)
    % Parameter vector is
    % theta =
    % (beta^1_0,            % Constant coefficient for firm 1
    % beta^1_1,             % Linear coefficient for firm 1
    % (beta^2_0,            % Constant coefficient for firm 2
    % beta^2_1,             % Linear coefficient for firm 2
    % delta^1_0,            % Constant competition effect for firm 1
    % delta^1_1,            % Linear competition effect for firm 1
    % delta^2_0,            % Constant competition effect for firm 2
    % delta^2_1)            % Linear competition effect for firm 2
    % rho                   % Correlation in Bivariate Normal Distribution

    % PRELIMINARY
    dim_p = KMSoptions.dim_p;
    suppX = KMSoptions.suppX;                   % Support for covariates [x1,x2], which is [-1 -1; -1 1 ;  1 -1 ; 1  1 ];
    psuppX = KMSoptions.psuppX;                 % Probability of support point occuring [0.1;0.2;0.3;0.4];
    dim_suppX = size(suppX,1);                  % Number of support points
    Dg_ineq = zeros(J1,dim_p);                     % Preset moment inequalities and
    Dg_eq = zeros(J2,dim_p);                       % moment equalities

    % Extract parameters (easier to read)
    beta1 = theta(1:2);
    beta2 = theta(3:4);
    delta1 = theta(5:6);
    delta2 = theta(7:8);
    rho    = theta(9);
    mu   = [0,0];
    sigma_rho = [1,rho;rho,1];
    r= -rho;
    sigma_r = [1,r; r,1]; 
        
    % For each point in the support, we get the theoretical gradient for each
    % potential point in the support
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

        % Probability of this support occuring
        pX= psuppX(ii);
        
        %%%% Gradient WRT beta1,beta2,delta1,delta2 %%%%%%%%%%%%%%%%%%%%%%%
        % Gradients for moments inequalities for entry game (See Pg 34, eq 5.3)
        % NB: g3 has r correlation.
        Dg3b1c = x1.*normpdf(-x1*(beta1+delta1)) .*normcdf((x2*beta2-r*(-x1*(beta1+delta1)))/sqrt(1-r^2))*pX;
        Dg3b2c = -x2.*normpdf(x2*beta2).*normcdf((-x1*(beta1+delta1)-r*(x2*beta2))/sqrt(1-r^2))*pX;
        Dg3d1c = x1.*normpdf(-x1*(beta1+delta1)) .*normcdf((x2*beta2-r*(-x1*(beta1+delta1)))/sqrt(1-r^2))*pX;
        Dg3d2c = zeros(1,2);
        
        % Gradients for moments inequalities for entry game (See Pg 34, eq 5.4)
        % NB: g4 has three copulas.  The first is r correlated.  The second two are rho correlated.
        Dg4b1c = -x1.*normpdf(-x1*(beta1+delta1)).*normcdf((x2*(beta2+delta2)-r*(-x1*(beta1+delta1)))/sqrt(1-r^2))*pX ...
            +((-x1.*normpdf(-x1*beta1).*normcdf((-x2*(beta2+delta2)-rho*(-x1*beta1))/sqrt(1-rho^2)))...
            -(-x1.*normpdf(-x1*beta1).*normcdf((-x2*beta2-rho*(-x1*beta1))/sqrt(1-rho^2))))*pX;
        Dg4b2c = x2.*normpdf(x2*(beta2+delta2)).*normcdf((-x1*(beta1+delta1)-r*(x2*(beta2+delta2)))/sqrt(1-r^2))*pX ...
            +((-x2.*normpdf(-x2*(beta2+delta2)).*normcdf((-x1*beta1-rho*(-x2*(beta2+delta2)))/sqrt(1-rho^2)))...
            -(-x2.*normpdf(-x2*beta2).*normcdf((-x1*beta1-rho*(-x2*beta2))/sqrt(1-rho^2))))*pX;
        Dg4d1c = -x1.*normpdf(-x1*(beta1+delta1)).*normcdf((x2*(beta2+delta2)-r*(-x1*(beta1+delta1)))/sqrt(1-r^2))*pX;
        Dg4d2c = x2.*normpdf(x2*(beta2+delta2)).*normcdf((-x1*(beta1+delta1)-r*(x2*(beta2+delta2)))/sqrt(1-r^2))*pX ...
            -x2.*normpdf(-x2*(beta2+delta2)).*normcdf((-x1*beta1-rho*(-x2*(beta2+delta2)))/sqrt(1-rho^2))*pX;
    
        % Gradients for moment equalities for entry game (See Pg 34, eq 5.1)
        Dg1b1c = x1.*normpdf(-x1*beta1).*normcdf((-x2*beta2-rho*(-x1*beta1))/sqrt(1-rho^2))*pX;
        Dg1b2c = x2.*normpdf(-x2*beta2).*normcdf((-x1*beta1-rho*(-x2*beta2))/sqrt(1-rho^2))*pX;
        Dg1d1c = zeros(1,2);
        Dg1d2c = zeros(1,2);

        % Gradients for moment equalities for entry game (See Pg 34, eq 5.2)
        Dg2b1c = -x1.*normpdf(x1*(beta1+delta1)).*(normcdf((x2*(beta2+delta2)-rho*(x1*(beta1+delta1)))/sqrt(1-rho^2)))*pX;
        Dg2b2c = -x2.*normpdf(x2*(beta2+delta2)).*(normcdf((x1*(beta1+delta1)-rho*(x2*(beta2+delta2)))/sqrt(1-rho^2)))*pX;
        Dg2d1c = -x1.*normpdf(x1*(beta1+delta1)).*(normcdf((x2*(beta2+delta2)-rho*(x1*(beta1+delta1)))/sqrt(1-rho^2)))*pX;
        Dg2d2c = -x2.*normpdf(x2*(beta2+delta2)).*(normcdf((x1*(beta1+delta1)-rho*(x2*(beta2+delta2)))/sqrt(1-rho^2)))*pX;

        %%%% Gradient WRT rho %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gradients for moments inequalities w.r.t. rho
        sgn_r   = -1;
        upper3 = [-x1*(beta1+delta1),x2*beta2];
        Dg3rho = -sgn_r*get_drho(r,mu,sigma_r,upper3)*pX;
        
        % Gradients for moments inequalities w.r.t. rho
        upper4_1 = [-x1*(beta1+delta1),x2*(beta2+delta2)];
        upper4_2 = [-x1*beta1,-x2*(beta2+delta2)];
        upper4_3 = [-x1*beta1,-x2*beta2];
        Dg4rho = sgn_r*get_drho(r,mu,sigma_r,upper4_1)*pX...
            +  get_drho(rho,mu,sigma_rho,upper4_2)*pX...
            -  get_drho(rho,mu,sigma_rho,upper4_3)*pX; 
        
        upper1 = [-x1*beta1,-x2*beta2];
        Dg1rho = -get_drho(rho,mu,sigma_rho,upper1)*pX;
        
        upper2 = [x1*(beta1+delta1),x2*(beta2+delta2)];
        Dg2rho = -get_drho(rho,mu,sigma_rho,upper2)*pX;

        %%%% Compile Moments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dg_ineq((ii-1)*2 + 1,:) = [Dg3b1c Dg3b2c Dg3d1c Dg3d2c Dg3rho];
        Dg_ineq((ii-1)*2 + 2,:) = [Dg4b1c Dg4b2c Dg4d1c Dg4d2c Dg4rho];
        Dg_eq((ii-1)*2 + 1,:)   = [Dg1b1c Dg1b2c Dg1d1c Dg1d2c Dg1rho];
        Dg_eq((ii-1)*2 + 2,:)   = [Dg2b1c Dg2b2c Dg2d1c Dg2d2c Dg2rho];

        % NOTE: These gradients are indeed J1 or 2*J2-by-dim_p, since G3b1 ect
        % is 2-by-1 (x1 and x2 include constant).
    end
    % Concatenate the positive and negative of Dmom_eq to transform into moment
    % inequalities
    Dg_eq = [Dg_eq ; -Dg_eq];

    % The way that the moment functions are written, m = f + g, g enters
    % additively.  The moment functions in this example have the incorrect
    % sign.  However, the gradient functions have the correct sign!
    
    % Sometimes get NaNs due to rounding precision.  Set equal to zero.
    Dg_ineq(isnan(Dg_ineq)) = 0; 
    Dg_eq(isnan(Dg_eq)) = 0; 
elseif DGP == 8   
    %% BCS Simulation
    % Example: BCS Entry Game
    % Page 20, Eq 5.1 in BCS

    % PRELIMINARY
    psuppX  = KMSoptions.psuppX;    % Prob. of X_q=1 for q=0,1,2,3.
    dX      = KMSoptions.dX;        % Number of covariates
    dim_p    = size(theta,1);       % Dimension of theta
    Dg_ineq = zeros(J1,dim_p);      % Preset moment inequalities and 
    Dg_eq = zeros(J2,dim_p);        % moment equalities

    % Extract parameters (easier to read)
    theta1 = theta(1);
    theta2 = theta(2);
    beta1 = theta(3);
    beta2 = theta(4);
    beta3 = theta(5);
    beta = [0; beta1; beta2; beta3];

    % GRADIENT COMPUTATION
    for qq = 1:dX
        px = psuppX(qq);

        % Gradient WRT theta1:
        Dg_ineq((qq-1)*2 + 1,1) = -(theta2-beta(qq))*px;                    % Moments m_{2q} in Eq 5.1 BCS
        Dg_ineq((qq-1)*2 + 2,1) = 0;                                        % Moments m_{3q} in Eq 5.1 BCS
        Dg_eq(qq,1)             = (1-theta2+beta(qq))*px;                   % Moments m_{1q} in Eq 5.1 BCS

        % Gradient WRT theta2:
        Dg_ineq((qq-1)*2 + 1,2) = (1-theta1+beta(qq))*px;                  % Moments m_{2q} in Eq 5.1 BCS
        Dg_ineq((qq-1)*2 + 2,2) = -px;                                     % Moments m_{3q} in Eq 5.1 BCS
        Dg_eq(qq,2)             = (1-theta1+beta(qq))*px;                  % Moments m_{1q} in Eq 5.1 BCS

        % Gradient WRT beta(qq), qq = 1,...,dX
        % beta(0) is normalized to 0, so its gradient is zero.  This is not
        % included in Dg_ineq and Dg_eq.
        if qq >= 2  
            % Gradient WRT beta(kk) is zero if qq ~= kk, since moment m_{i,qq}
            % depends only on beta(qq)
            Dg_ineq((qq-1)*2 + 1,1+qq) = (theta1+theta2-1-2*beta(qq))*px;          % Moments m_{2q} in Eq 5.1 BCS
            Dg_ineq((qq-1)*2 + 2,1+qq) = px;                                       % Moments m_{3q} in Eq 5.1 BCS
            Dg_eq(qq,1+qq)             = -(2 - theta1 - theta2 + 2*beta(qq))*px;    % Moments m_{1q} in Eq 5.1 BCS
        end
    end

    % Concatenate the positive and negative of Dmom_eq to transform into moment
    % inequalities
    Dg_eq = [Dg_eq ; -Dg_eq];
end
end

% Below are functions used for computing the gradient in DGP=7.
function drho = get_drho(rho,mu,sigma_rho,upper)
c1 = mexBVNcdf(upper,mu,sigma_rho)./((1-rho.^2).^2);
c2 = rho - rho.^3;
c3 = (1+rho.^2).*get_EXiXj(1,2,mu,sigma_rho,upper);
c4 = -rho.*(get_EXiXj(1,1,mu,sigma_rho,upper)+get_EXiXj(2,2,mu,sigma_rho,upper));
drho = c1.*(c2+c3+c4);
end

function EXiXj = get_EXiXj(i,j,mu,sigma_rho,upper)
if i==1 && j==2
    EXiXj = sigma_rho(1,2) ...
        + sigma_rho(2,1).*(-upper(:,1).*get_Fk(upper(:,1),1,mu,sigma_rho,upper))...
        + sigma_rho(1,2).*(-upper(:,2).*get_Fk(upper(:,2),2,mu,sigma_rho,upper))...
        + (sigma_rho(1,1).*sigma_rho(2,2)-sigma_rho(1,2).*sigma_rho(2,1))...
        .*get_Fqr(upper,mu,sigma_rho,upper);
elseif i==1 && j==1
    EXiXj= sigma_rho(1,1)...
        + sigma_rho(1,1).*(-upper(:,1).*get_Fk(upper(:,1),1,mu,sigma_rho,upper))...
        + sigma_rho(1,2).^2./sigma_rho(2,2)...
        .*(-upper(:,2).*get_Fk(upper(:,2),2,mu,sigma_rho,upper))...
        + (sigma_rho(1,2).*sigma_rho(1,1)-sigma_rho(1,2).^2.*sigma_rho(2,1)./sigma_rho(2,2))...
        .*get_Fqr(upper,mu,sigma_rho,upper);
elseif i==2 && j==2
    EXiXj= sigma_rho(2,2)...
        + sigma_rho(2,2).*(-upper(:,2).*get_Fk(upper(:,2),2,mu,sigma_rho,upper))...
        + sigma_rho(2,1).^2./sigma_rho(1,1)...
        .*(-upper(:,1).*get_Fk(upper(:,1),1,mu,sigma_rho,upper))...
        + (sigma_rho(2,1).*sigma_rho(2,2)-sigma_rho(2,1).^2.*sigma_rho(1,2)./sigma_rho(1,1))...
        .*get_Fqr(upper,mu,sigma_rho,upper);
end
end

function Fqr = get_Fqr(X,mu,sigma_rho,upper)
% X: n-by-2
% upper: n-by-2
% mean: 1-by-2
% sigma: 2-by-2
alpha = mexBVNcdf(upper,mu,sigma_rho);
Fqr = mvnpdf(X,mu,sigma_rho)./alpha;
end

function Fk = get_Fk(xn,i,mu,sigma_rho,upper)
n = length(xn);
C = sigma_rho;
A = inv(sigma_rho);
if i==1
    j=2;
elseif i==2
    j=1;
end
A_1 = A(j,j);
A_1_inv = inv(A_1);
C_1 = C(j,j);
c_nn = C(i,i);
c = C(j,i);
mu_1 = mu(j);
mu_n = mu(i);
f_xn = zeros(n,1);
p = mexBVNcdf(upper,mu,sigma_rho);
for l=1:n
    m = mu_1 + (xn(l) - mu_n) .* c/c_nn;
    f_xn(l) = exp(-0.5 .* (xn(l) - mu_n).^2./c_nn) .* normcdf(upper(l,j),m,sqrt(A_1_inv));
end
Fk = 1./p .* 1./sqrt(2 * pi * c_nn) .* f_xn;
end
