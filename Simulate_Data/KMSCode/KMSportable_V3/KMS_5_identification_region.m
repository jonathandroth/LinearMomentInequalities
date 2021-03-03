function [Identification_region] =  KMS_5_identification_region(theta_true,theta_0,LB_theta,UB_theta,A_theta,b_theta,KMSoptions);
%% Code Description:
% This function computes the true identification region for the games
% and the BCS example.

%% Extract information
DGP = KMSoptions.DGP;
dim_p = size(theta_true,1);
KMSoptions.dim_p = dim_p;
Identification_region = zeros(dim_p,2);
options_fmincon = KMSoptions.options_fmincon;
options_multistart  = KMSoptions.options_multistart;

%% Population moments
if DGP == 5 || DGP == 6 
    % True population moment:
    % PRELIMINARY    
    suppX = KMSoptions.suppX;                   % Support for covariates [x1,x2] = [-1 -1; -1 1 ;  1 -1 ; 1  1 ]; in baseline spec
    dim_suppX = size(suppX,1);                  % Number of support points
    J1 = 2*dim_suppX;
    J2 = 2*dim_suppX;
    selp = KMSoptions.selp;                     % Selection probability
    psuppX = KMSoptions.psuppX;
    KMSoptions.J1 = J1;
    KMSoptions.J2 = J2;
    
    % Extract parameters (easier to read)
    beta1 = theta_true(1:2);
    beta2 = theta_true(3:4);
    delta1 = theta_true(5:6);
    delta2 = theta_true(7:8);
    
    % Preset output variables
    f_ineq_pop = zeros(J1,1);                        % Preset moment inequalities and
    f_eq_pop = zeros(J2,1);                          % moment equalities

    % MOMENT COMPUTATION
    % The moments f_j(Y,X) is equal to the frequency of observing (Y,X):
    % (Y1=y1,Y2=y2,X1=x1,X2=x2).
    % Potential outcomes are:
    % (y1,y2) = (0,0)  (both firms do not enter)
    % (y1,y2) = (1,1)  (both firms enter)
    % (y1,y2) = (1,0)  (only firm 1 enters)
    % (y1,y2) = (0,1)  (Only firm 2 enters)
    for ii = 1:dim_suppX
        pX= psuppX(ii);

        % Pick support point (include constant)
        x1 = [1,suppX(ii,2)];
        x2 = [1,suppX(ii,4)];
        
        % Moment inequalities for entry game (See Pg 34, eq 5.3)
        f_ineq_pop((ii-1)*2 + 1,1) = normcdf(-x1*(beta1+delta1))*(1-normcdf(-x2*beta2))*pX;

        % Moment inequalities for entry game (See Pg 34, eq 5.4)
        f_ineq_pop((ii-1)*2 + 2,1) = ...
            -normcdf(-x1*(beta1+delta1))*(1-normcdf(-x2*(beta2+delta2)))*pX...
            -normcdf(-x1*beta1)*(normcdf(-x2*(beta2+delta2))-normcdf(-x2*beta2))*pX;
       
        % Moment equalities for entry game (See Pg 34, eq 5.1)
        f_eq_pop((ii-1)*2 + 1,1)   = normcdf(-x1*beta1)*normcdf(-x2*beta2)*pX;

        % Moment equalities for entry game (See Pg 34, eq 5.2)
        f_eq_pop((ii-1)*2 + 2,1) =(1-normcdf(-x1*(beta1+delta1)))*(1-normcdf(-x2*(beta2+delta2)))*pX;   
        
        % Adjust moment inequalities using selp to get moment equalities
        width = f_ineq_pop((ii-1)*2 + 1,1) + f_ineq_pop((ii-1)*2 + 2,1);
        f_ineq_pop((ii-1)*2 + 1,1) = f_ineq_pop((ii-1)*2 + 1,1) - width*selp;
        f_ineq_pop((ii-1)*2 + 2,1) = -f_ineq_pop((ii-1)*2 + 1,1);
    end
    % Concatenate the positive and negative of f_eq to transform into moment
    % inequalities
    f_eq_pop =[f_eq_pop;-f_eq_pop];
elseif DGP == 7
    % True population moment:
    % PRELIMINARY    
    suppX = KMSoptions.suppX;                   % Support for covariates [x1,x2] = [-1 -1; -1 1 ;  1 -1 ; 1  1 ]; in baseline spec
    dim_suppX = size(suppX,1);                  % Number of support points
    J1 = 2*dim_suppX;
    J2 = 2*dim_suppX;
    selp = KMSoptions.selp;                     % Selection probability
    psuppX = KMSoptions.psuppX;
    KMSoptions.J1 = J1;
    KMSoptions.J2 = J2;
    
    % Extract parameters (easier to read)
    beta1 = theta_true(1:2);
    beta2 = theta_true(3:4);
    delta1 = theta_true(5:6);
    delta2 = theta_true(7:8);
    rho     = theta_true(9);
    
    % Preset output variables
    f_ineq_pop = zeros(J1,1);                        % Preset moment inequalities and
    f_eq_pop = zeros(J2,1);                          % moment equalities

    % MOMENT COMPUTATION
    % The moments f_j(Y,X) is equal to the frequency of observing (Y,X):
    % (Y1=y1,Y2=y2,X1=x1,X2=x2).
    % Potential outcomes are:
    % (y1,y2) = (0,0)  (both firms do not enter)
    % (y1,y2) = (1,1)  (both firms enter)
    % (y1,y2) = (1,0)  (only firm 1 enters)
    % (y1,y2) = (0,1)  (Only firm 2 enters)
    
 
    for ii = 1:dim_suppX
        pX= psuppX(ii);

        % Pick support point (include constant)
        x1 = [1,suppX(ii,2)];
        x2 = [1,suppX(ii,4)];
        
        % Moment inequalities for entry game (See Pg 34, eq 5.3)
        f_ineq_pop((ii-1)*2 + 1,1) = mexBVNcdf([-x1*(beta1+delta1),x2*beta2],[0;0],[1,-rho;-rho,1])*pX;

        % Moment inequalities for entry game (See Pg 34, eq 5.4)
        f_ineq_pop((ii-1)*2 + 2,1) = ...
        -mexBVNcdf([-x1*(beta1+delta1),x2*(beta2+delta2)],[0;0],[1,-rho;-rho,1])*pX...
        -mexBVNcdf([-x1*beta1,-x2*(beta2+delta2)],[0;0],[1,rho;rho,1])*pX ...
        +mexBVNcdf([-x1*beta1,-x2*beta2],[0;0],[1,rho;rho,1])*pX;
           
        % Moment equalities for entry game (See Pg 34, eq 5.1)
        f_eq_pop((ii-1)*2 + 1,1)   =  mexBVNcdf([-x1*beta1,-x2*beta2],[0;0],[1,rho;rho,1])*pX;

        % Moment equalities for entry game (See Pg 34, eq 5.2)
        f_eq_pop((ii-1)*2 + 2,1) =mexBVNcdf([x1*(beta1+delta1),x2*(beta2+delta2)],[0;0],[1,rho;rho,1])*pX;

        % Adjust moment inequalities using selp to get moment equalities
        width = f_ineq_pop((ii-1)*2 + 1,1) + f_ineq_pop((ii-1)*2 + 2,1);
        f_ineq_pop((ii-1)*2 + 1,1) = f_ineq_pop((ii-1)*2 + 1,1) - width*selp;
        f_ineq_pop((ii-1)*2 + 2,1) = -f_ineq_pop((ii-1)*2 + 1,1);
    end
    % Concatenate the positive and negative of f_eq to transform into moment
    % inequalities
    f_eq_pop =[f_eq_pop;-f_eq_pop];
    
elseif DGP == 8  
    % True population moment:
    % PRELIMINARY
    dX = KMSoptions.dX;
    psuppX = KMSoptions.psuppX;
    selp = KMSoptions.selp;
    J1 = 2*dX;
    J2 = dX;
    J = J1 + 2*J2;
    KMSoptions.J1 = J1;
    KMSoptions.J2 = J2;
    KMSoptions.J = J;

    % Preset output variables
    f_ineq_pop = zeros(J1,1);                        % Preset moment inequalities and
    f_eq_pop = zeros(J2,1);                         % moment equalities

    % Extract parameters (easier to read)
    theta1_true = theta_true(1);
    theta2_true = theta_true(2);
    beta1_true = theta_true(3);
    beta2_true = theta_true(4);
    beta3_true = theta_true(5);
    beta_true = [0; beta1_true; beta2_true; beta3_true];

    % MOMENT COMPUTATION
    for ii = 1:dX
        f_ineq_pop((ii-1)*2 + 2,1) = (1-theta1_true+beta_true(ii))*(theta2_true-beta_true(ii)) ...
                                    + (1-selp)*(theta1_true-beta_true(ii))*(theta2_true-beta_true(ii)); % Moments m_{3q} in Eq 5.1 BCS
        f_ineq_pop((ii-1)*2 + 1,1) = -f_ineq_pop((ii-1)*2 + 2,1);                                   % Moments m_{2q} in Eq 5.1 BCS
        f_eq_pop(ii,1)             = (1-theta1_true+beta_true(ii))*(1-theta2_true+beta_true(ii));   % Moments m_{1q} in Eq 5.1 BCS

        
        
        
        % Scale by psuppX(ii) to get correct freqency
        f_ineq_pop((ii-1)*2 + 1,1) = f_ineq_pop((ii-1)*2 + 1,1)*psuppX(ii);
        f_ineq_pop((ii-1)*2 + 2,1) = f_ineq_pop((ii-1)*2 + 2,1)*psuppX(ii);
        f_eq_pop(ii,1)             = f_eq_pop(ii,1)*psuppX(ii);
    end

    % Concatenate the positive and negative of f_eq to transform into moment
    % inequalities
    f_eq_pop = [f_eq_pop ; -f_eq_pop];
end


%% Identification region
% We solve min/max p'theta in all basis vector directions
options_fmincon.TolCon = 1e-10;
options_fmincon.TolFun = 1e-10;
options_fmincon.TolX = 1e-10;

options_multistart.Display = 'off';
mult_num = 50;

for ii = 1:dim_p
   fprintf('Starting Identification Region Component=%d \n',ii)

   % maximize direction
   p = zeros(dim_p,1);
   p(ii) = 1;
   
   % Objective function and constraint:
   obj_true =  @(theta)KMS_51_IRobjective(theta,p);
   constraint_true =  @(theta)KMS_52_IRconstraint(theta,f_ineq_pop,f_eq_pop,KMSoptions);
    
   % Solve using fmincon from each initial point theta_0_fminimax. 
   %  NB: we select the scale option in fmincon to avoid scaling issues.
  % [x,fval,exitflag] = fmincon(obj_true,theta_true,A_theta,b_theta,[],[],...
  %                  LB_theta,UB_theta,constraint_true,options_fmincon);

   problem = createOptimProblem('fmincon','x0',theta_0,...
        'objective',obj_true,'Aineq', A_theta, 'bineq', b_theta,'lb',LB_theta,'ub',UB_theta,...
        'nonlcon',constraint_true,'options',options_fmincon);
       
    [~,fval,exitflag] = run(options_multistart,problem,mult_num);
    if exitflag > 0
        Identification_region(ii,2) = -fval;
    else
        Identification_region(ii,2) = nan;
    end

   % minimize direction
   p = zeros(dim_p,1);
   p(ii) = -1;
   
   % Objective function and constraint:
   obj_true =  @(theta)KMS_51_IRobjective(theta,p);
   constraint_true =  @(theta)KMS_52_IRconstraint(theta,f_ineq_pop,f_eq_pop,KMSoptions);
   
   problem = createOptimProblem('fmincon','x0',theta_0,...
        'objective',obj_true,'Aineq', A_theta, 'bineq', b_theta,'lb',LB_theta,'ub',UB_theta,...
        'nonlcon',constraint_true,'options',options_fmincon);

    [~,fval,exitflag] = run(options_multistart,problem,mult_num);
    if exitflag > 0
        Identification_region(ii,1) = fval;
    else
        Identification_region(ii,1) = nan;
    end
end



end