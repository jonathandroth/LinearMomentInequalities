function [theta_hat,theta_optbound,c,CV,EI,flag_opt] =  KMS_3_EAM(q,sgn_q,theta_feas,theta_init,c_init,CV_init,maxviol_init,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions)
%% Code description: EAM
%  This function executes the EAM algorithm and outputs the KMS interval. 
%  For a complete description of EAM see Pages 12-13.
%
%  This function outputs the value that solves problem Eq 2.16, Pg 13:
%
%  (Eq 2.16)
%  max/min_{p} p'theta
%  s.t. sqrt(n)(moments(theta))/sigma(theta) <= c(theta)
%
%  If sgn_q = 1, then we solve max_p p'theta. If sgn_q = -1, then we solve
%  min_p p'theta.
%
% INPUTS:
%   q                   dim_p-by-1 directional vector.  This is either
%                       p or -p
%
%   sgn_q               Equal to 1 or -1.  This determines the value of the
%                       problem.  Since we are solving min/max p'theta, the
%                       value of the problem is sgn_q*q.'*theta_opt
%
%   theta_feas          Feasible point found in auxiliary search
%
%   f_ineq,f_eq         Empirical moments
%
%   f_ineq_keep,f_eq_keep      Moments to keep   
%
%  f_stdev_ineq,f_stdev_eq     Standard deviation of empirical moments
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
%   theta_hat           dim_p-by-1 parameter vector that solves Eq 2.16
%
%   theta_optbound      1-by-1 optimal value/projection: abs(q')*theta_hat.
%   
%   flag_opt            flag_opt = 1 if program converged
%                       flag_opt = -1 if no feasible points can be found
%                       (output point that minimizes constraint violation)
%                       flag_opt = 0 if EAM did not converge in the maximum
%                       number of iterations.
%   c                   Critical value c(theta)
%
%   CV                  Constriant violation -- 0 if converged
%
%   EI                  Expected improvement -- should be small if
%                       converged

%% Extract relevant information from KMSoptions
LB_theta            = KMSoptions.LB_theta;
UB_theta            = KMSoptions.UB_theta;
A_theta             = KMSoptions.A_theta;
b_theta             = KMSoptions.b_theta;
CI_lo               = KMSoptions.CI_lo;
CI_hi               = KMSoptions.CI_hi;
sample_method       = KMSoptions.sample_method;
dim_p               = KMSoptions.dim_p;
dace_theta          = KMSoptions.dace_theta;
dace_lob            = KMSoptions.dace_lob;
dace_upb            = KMSoptions.dace_upb;
e_points_init       = KMSoptions.e_points_init;
options_fmincon     =  KMSoptions.options_fmincon;
options_linprog     = KMSoptions.options_linprog;
EAM_maxit           = KMSoptions.EAM_maxit;
EAM_minit           = KMSoptions.EAM_minit;
EAM_tol             = KMSoptions.EAM_tol;
h_rate              = KMSoptions.h_rate;
h_rate2             = KMSoptions.h_rate2;
parallel            = KMSoptions.parallel;
unif_num            = KMSoptions.unif_num;
EI_num              = KMSoptions.EI_num;
EI_multi_num        = KMSoptions.EI_multi_num;
EAM_obj_tol         = KMSoptions.EAM_obj_tol;
r_min               = KMSoptions.r_min;
EAM_maxviol_tol     = KMSoptions.EAM_maxviol_tol;
EAM_thetadistort    = KMSoptions.EAM_thetadistort;

%% Initialize EAM
% Number of feasible points
num_feas = size(theta_feas,1);

% Draw initial points uniformly from Theta 
if sample_method == 0 || sample_method == 1 
    % Bounds are a box, so use Latin hypercute sampling
    theta_add = KMS_AUX2_drawpoints(e_points_init,dim_p,LB_theta,UB_theta,KMSoptions);
else
    % Bounds define a polytope, so use hit-and-run sampling
    theta_add = KMS_AUX2_drawpoints(e_points_init,dim_p,LB_theta,UB_theta,KMSoptions,A_theta,b_theta,theta_feas.');
end

% Set theta_Estep (initial set of points for EAM) to be the union of the
% uniformly drawn points, bounds on Theta, and theta_feas.  
theta_Estep = [theta_add ; theta_feas];

% theta_Astep is the set of thetas that EAM has explored
% So that the program runs correctly, I initiate it to be the empty set.
% theta_Astep will be updated shortly in the EAM algorithm
theta_Astep     = theta_init;
c_Astep         = c_init;
CV_Astep        = CV_init;
maxviol_Astep   = maxviol_init;

% Change in optimal projection is undefined for the first iteration
opt_val_old = nan;          

% Contraction counter counts the degree to which we contract the parameter
% space.  It is an integer with value {0,1,2,...,EAM_maxit}.
% If on the previous iteration we failed to find a feasible point, then we
% contract the parameter space by increasing the counter.
% If, on the other hand, the new feasible point is close to the boundary,
% we expand the parameter space by reducing the counter.
% The parameter space is contracted/expanded at rate 1/h_rate^{counter}.
contraction_counter = 0;

% opt_bound is upper bound on the parameter space in direction q.
% opt_dagger is the upper bound on the contracted parameter space in
% direction q.  opt_dagger changes iteration-to-iteration.
if sgn_q == 1
    opt_dagger   = CI_hi;
    opt_bound    = CI_hi;
else
    opt_dagger   = CI_lo;
    opt_bound    = CI_lo;
end

%% EAM optimization routine
% We run EAM for up to EAM_maxit times.  Each loop adds more evaluation points
% until we yield convergence.
fprintf('Iteration     |     Opt Proj     | Change in EI proj        |     Change     |     Max Violation   | Feasible points    |  Multi. Start Num.  | Percent conv. EI>0  |  Contraction counter \n')
fprintf('------------------------------------------------------------------------------------------------------------------------------------------------------- \n')
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
    
    % Step 2) A-Step 
    % This step interpolates critical values outside of the evaluation 
    % points
    % Make sure design points are not too close together
    [theta_dmodel,ind] = uniquetol(theta_Astep,1e-10,'ByRows',true); 
    c_dmodel = c_Astep(ind,:);
    dmodel = dacefit(theta_dmodel,c_dmodel,@regpoly0,@corrgauss,dace_theta,dace_lob,dace_upb);
    
    % Step 3) M-Step 
    % This step draws new point for next iteration
    % The next point(s) are drawn using Jones' expected improvement method
    % with constraints.  Briefly: a prior is put over the constraints and
    % the next point(s) are drawn to maximize the expected gain in the
    % objective function.  The expected improvement function can be written
    % as a minimax problem.  We, however, write the minimax problem so that 
    % fmincon can solve it, since we found that fminimax is unstable in
    % simulations.
    %
    % Define h_j(theta) = sqrt(n)*m_j(X,theta)/sigma(X), which are the
    % the standardized moments.  The expected improvement objective 
    % function can be written as:
    %   
    %   max_j(q'theta - q'theta_#)_{+}*Phi( (h_j(theta)  -
    %   c_L(theta))/sigma_L(theta))
    %
    % where  c_L(theta) and sigma_L(theta) are estimated from the DACE
    % model, and theta_# is the point that maximizes the objective function 
    % q'theta subject to the constraint that theta is in the set
    % S = {theta_1,...,theta_L : CV(theta_l) = 0},
    % i.e., the set of theta's already explored that are feasible.
    % We can further simplify the problem by searching over
    % the space of Theta satisfying q'theta  >= q'theta_# and drop the
    % (.)_{+} operator.  
    % Find points that have 0 constraint violation.  The auxiliary
    % feasible search gaurantees that such a point exists (namely,
    % theta_feas.
    feas = find(CV_Astep == 0);
    theta_feas = theta_Astep(feas,:);
    maxviol_feas = maxviol_Astep(feas);
    [~,ind] = max(theta_feas*q);
    theta_hash = theta_feas(ind,:).';
    maxviol_hash = maxviol_feas(ind);
    
    % Linear constraints:
    % We require that q'theta >= q'theta_#
    % Constraints are in the form D*theta <= d.  So set d=-q'theta_# and
    %   D = -q.'.
    % If the projection vector, q, is a basis vector we can embed these
    % constraints into box constraints.  Otherwise, we embed them as a
    % polytope constraints Ax <= b.
    if sample_method == 0
        % Update lower/upper bounds
        LB_EI = LB_theta;
        UB_EI = UB_theta;
        if sgn_q == 1
            LB_EI(KMSoptions.component) = q.'*theta_hash;
            UB_EI(KMSoptions.component) = q.'*theta_hash + (opt_bound-q.'*theta_hash)/(h_rate^contraction_counter); % Shrink parameter space by r_tate
        else
            LB_EI(KMSoptions.component) = abs(q).'*theta_hash - (abs(q).'*theta_hash - opt_bound)/(h_rate^contraction_counter);
            UB_EI(KMSoptions.component) = abs(q).'*theta_hash;
        end
        
        % Due to numerical error, it is possible that the LB or UB violates
        % the UB and LB imposed by user.  We correct for this by
        % overwriting LB_EI and UB_EI if either violates LB or UB.
        LB_EI(KMSoptions.component) = max(CI_lo,LB_EI(KMSoptions.component));
        UB_EI(KMSoptions.component) = min(CI_hi,UB_EI(KMSoptions.component));
        
        % In this case there are no polytope constraints Ax <= b.
        A_EI = [];
        b_EI = [];
    elseif sample_method == 1 || sample_method == 2
         % Update lower/upper bounds
        LB_EI = LB_theta;
        UB_EI = UB_theta;
        
        % Update lower/upper bounds Ax <= b.
        A_EI = A_theta;
        b_EI = b_theta;
        
        % Include bound that q.'theta_hash <= q.'theta.
        A_EI = [A_EI ; -q.'];  
        b_EI = [b_EI ; -q.'*theta_hash];
        
        % Include constraint q.'theta <= q.'theta_{bound}*scale
        % where theta_{bound} is the upper bound on the parameter space in 
        % direction q, and scale is a constant that depends on the
        % contraction counter and the distance between the upper bound and
        % theta_hash.
        d = abs(opt_bound - q.'*theta_hash);
        A_EI = [A_EI ; q.'];
        b_EI = [b_EI ; q.'*theta_hash + d/(h_rate^contraction_counter)];
    elseif sample_method == 3
        % Update lower/upper bounds
        LB_EI = LB_theta;
        UB_EI = UB_theta;
        if sgn_q == 1
            LB_EI(KMSoptions.component) = q.'*theta_hash;
            UB_EI(KMSoptions.component) = q.'*theta_hash + (opt_bound-q.'*theta_hash)/(h_rate^contraction_counter); % Shrink parameter space by r_tate
        else
            LB_EI(KMSoptions.component) = abs(q).'*theta_hash - (abs(q).'*theta_hash - opt_bound)/(h_rate^contraction_counter);
            UB_EI(KMSoptions.component) = abs(q).'*theta_hash;
        end
        
        % Due to numerical error, it is possible that the LB or UB violates
        % the UB and LB imposed by user.  We correct for this by
        % overwriting LB_EI and UB_EI if either violates LB or UB.
        LB_EI(KMSoptions.component) = max(CI_lo,LB_EI(KMSoptions.component));
        UB_EI(KMSoptions.component) = min(CI_hi,UB_EI(KMSoptions.component));
        
        % Bound transform if draw-and-discard method is used
        [LB_EI,UB_EI] = bound_transform(LB_EI,UB_EI,KMSoptions);

        % Update lower/upper bounds Ax <= b.
        A_EI = A_theta;
        b_EI = b_theta;     
    end
    
    % Update opt_dagger
    % We have contracted the parameter space, so we need to update the
    % maximum value of q'theta s.t. theta in parameter space.
    [theta_dagger,opt_dagger] =  linprog(-q,A_EI,b_EI,[],[],LB_EI,UB_EI,[],options_linprog);
    opt_dagger = -opt_dagger;

    % Draw initial points between evaluation points.
    % The points are drawn in a particular way so that the search algorithm
    % is more likely to converge.
    % r_max is the distance from theta_hash to the boundary.
    r_max = abs(opt_dagger - q.'*theta_hash);
    r_max = max(r_max, r_min);
    [theta_keep,EI_keep]  = KMS_36_drawpoints(theta_hash,q,r_max,r_min,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,dmodel,LB_EI,UB_EI,A_EI,b_EI,KMSoptions);
    
    if ~isempty(theta_keep)
        % Draw initial points between evaluation points.
        % The points are drawn in a particular way so that the search algorithm
        % is more likely to converge.
        % (See Matthias Schonlau; William J Welch; Donald R Jones, 1998)
        theta_0_fminimax = KMS_AUX4_MSpoints([theta_keep;theta_hash.']);
        theta_0_fminimax = unique(theta_0_fminimax,'rows');

        % Find thetas that have positive EI and drop those with EI = 0
        Eimprovement = @(theta)KMS_37_EI_value(theta,q,theta_hash,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,dmodel,KMSoptions);
        EI_fminimax = zeros(size(theta_0_fminimax,1),1);
        if parallel
            parfor jj = 1:size(theta_0_fminimax,1)
                EI_fminimax(jj,1) = -(max(Eimprovement( theta_0_fminimax(jj,:).')));
            end
        else
            for jj = 1:size(theta_0_fminimax,1)
                EI_fminimax(jj,1) = -(max(Eimprovement( theta_0_fminimax(jj,:).')));
            end
        end
        % Keep solutions with positive expected improvement
        ind = find(EI_fminimax <= 0);
        theta_0_fminimax(ind,:) = [];
        EI_fminimax(ind,:) = [];

        % Sort by EI
        [EI_fminimax,I] = sort(EI_fminimax,'descend');
        theta_0_fminimax = theta_0_fminimax(I,:);

        % Keep top EI_multi_num
        theta_0_fminimax(EI_multi_num+1:end,:) = [];
    else
        theta_0_fminimax = [];
    end
    % Include theta_hash and theta_eps
    theta_eps = theta_hash.' + q.'*(1e-4);   
    theta_0_fminimax = [theta_0_fminimax;theta_hash.';theta_eps];
    
    % Number of initial points:
    multistart_num = size(theta_0_fminimax,1);
    
    % Run fmincon with multistart
    theta_Mstep      = zeros(multistart_num,dim_p);
    EI_Mstep         = zeros(multistart_num,1);
    flag_conv        = zeros(multistart_num,1);
          
    % Objective and constraint
    objective_Eimprovement = @(theta)KMS_34_EI_objective(theta,KMSoptions);
    constraint_Eimprovement = @(theta)KMS_35_EI_constraint(theta,q,theta_hash,...
            f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,dmodel,KMSoptions);    
        
    % Solve using fmincon from each initial point theta_0_fminimax. 
    if parallel
        parfor ii = 1:multistart_num
            try
                theta_aug = [theta_0_fminimax(ii,:).';0];

                [x,fval,exitflag] = fmincon(objective_Eimprovement,theta_aug,[A_EI, zeros(size(A_EI,1),1)],b_EI,[],[],...
                                [LB_EI;-inf],[UB_EI;inf],constraint_Eimprovement,options_fmincon);
                theta_Mstep(ii,:) = x(1:dim_p,1).';
                EI_Mstep(ii,1) =  -fval;
                flag_conv(ii,1) = exitflag;
            catch
                theta_Mstep(ii,:) = theta_0_fminimax(ii,:);
                EI_Mstep(ii,1)    = 0;
                flag_conv(ii,1)   = -1;
            end
        end
    else
        for ii = 1:multistart_num
            try
                theta_aug = [theta_0_fminimax(ii,:).';0];

                [x,fval,exitflag] = fmincon(objective_Eimprovement,theta_aug,[A_EI, zeros(size(A_EI,1),1)],b_EI,[],[],...
                                [LB_EI;-inf],[UB_EI;inf],constraint_Eimprovement,options_fmincon);
                theta_Mstep(ii,:) = x(1:dim_p,1).';
                EI_Mstep(ii,1) =  -fval;
                flag_conv(ii,1) = exitflag;
            catch
                theta_Mstep(ii,:) = theta_0_fminimax(ii,:);
                EI_Mstep(ii,1)    = 0;
                flag_conv(ii,1)   = -1;
            end
        end
    end
    % Keep soltuions that are feasible 
    ind = find(flag_conv<= 0);
    theta_Mstep(ind,:) = [];
    EI_Mstep(ind,:) = [];
    
    % Check solutions are inside the parameter space
    A_aug = [A_EI ; eye(dim_p) ; -eye(dim_p)];
    b_aug = [b_EI ; UB_EI ; -LB_EI];
    size_opt = size(theta_Mstep,1);
    ind = find(max(A_aug*(theta_Mstep.') - repmat(b_aug,[1,size_opt])) > 0).';
    theta_Mstep(ind,:) = [];
    EI_Mstep(ind,:) = [];
    
    % Percent of runs that converged to theta with positive EI.
    percent_conv = 100*size(find(EI_Mstep>1e-15),1)/multistart_num;

    % Drop solutions with expected improvement = 0
    ind = find(EI_Mstep ~= 0);
    theta_Mstep = theta_Mstep(ind,:);
    EI_Mstep    = EI_Mstep(ind);

    % Sort by expected improvement.
    % NB: we include both initial points and those that were found from the
    % maximization problem
    EI_Mstep = [EI_Mstep;EI_keep];
    theta_Mstep = [theta_Mstep;theta_keep];
    
    [EI_Mstep,I] = sort(EI_Mstep,'descend');
    theta_Mstep = theta_Mstep(I,:);
  
    [theta_Mstep,I2] = uniquetol(theta_Mstep,1e-8,'ByRows',true); 
    EI_Mstep = EI_Mstep(I2);
    
    % Resort (problem with shuffling after uniquetol)
    [EI_Mstep,I] = sort(EI_Mstep,'descend');
    theta_Mstep = theta_Mstep(I,:);
    
    % Keep top EI_num points
    theta_Mstep(EI_num+1:end,:) = [];
    EI_Mstep(EI_num+1:end,:) = [];
    
   % Plus draw unif_num from {theta : p'theta >= p'theta#}
   if sample_method == 0
       theta_draw = KMS_AUX2_drawpoints(unif_num,dim_p,LB_EI,UB_EI,KMSoptions);
   else
       theta_draw = KMS_AUX2_drawpoints(unif_num,dim_p,LB_EI,UB_EI,KMSoptions,A_EI,b_EI,theta_hash);
   end

   theta_Estep = [theta_Mstep;theta_draw];
   if isempty(I) ==0
        EI = EI_Mstep(1);
   else
        EI = nan;  
   end
   
    % Also add a small distortion of theta# to theta_Estep
    delta1 = abs(maxviol_hash)/(h_rate2^contraction_counter);
    delta2 = EAM_thetadistort;
    if sgn_q == 1
        theta_eps1 = min(theta_hash.' + q.'*delta1,UB_theta.');
        theta_eps2 = min(theta_hash.' + q.'*delta2,UB_theta.');  
    else
        theta_eps1 = max(theta_hash.' + q.'*delta1,LB_theta.'); 
        theta_eps2 = max(theta_hash.' + q.'*delta2,LB_theta.');  
    end
    theta_eps = [theta_eps1;theta_eps2];
    
    % Check theta_eps1,theta_eps2 are inside the parameter space
    A_aug = [A_EI ; eye(dim_p) ; -eye(dim_p)];
    b_aug = [b_EI ; UB_EI ; -LB_EI];
    size_opt = size(theta_eps,1);
    ind = find(max(A_aug*(theta_eps.') - repmat(b_aug,[1,size_opt])) > 0).';
    theta_eps(ind,:) = [];
    theta_Estep = [theta_Estep;theta_eps];

    % Step 4) Print Results and Convergence
    % Program converges when expected improvement is less than
    % "best" current value of objective function divided by 100.
    opt_val = sgn_q*q.'*theta_hash;
    if isempty(I) ==0
        opt_EI_proj = sgn_q*q.'*(theta_Mstep(1,:).'); 
    else
        opt_EI_proj = nan;
    end
    change_EI_proj = abs(opt_EI_proj  - opt_val);
    change_proj =abs(opt_val - opt_val_old);
    feas_points = sum(CV_Astep == 0);
    Output = [iter, opt_val, change_EI_proj, change_proj,maxviol_hash,feas_points,multistart_num , percent_conv,contraction_counter];
    fprintf('%9.4f     |   %9.4f      | %9.4e               | %9.4f      | %9.4f    | %9.4f          | %9.4f          | %9.4f           | %9.4f            \n',Output)

    % Check for convergence
    % If the best feasible point are too close to the parameter boundary,
    % conclude that we have converged but output warning -- the KMS theory
    % does not hold if the parameter is on the boundary, so we may not get
    % correct covereage
    if abs(opt_val - opt_bound) < 1e-4
        theta_hat     = theta_hash;
        theta_optbound= opt_val;
        [c,CV] = KMS_31_Estep(theta_hash.',f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
        EI =  EI(1);
        flag_opt =1;
        warning('Parameter is on the boundary.  The confidence set might not deliver the correct coverage.  Consider expanding the parameter space.')
        return;
    end
    if (iter >= EAM_minit &&  change_EI_proj < EAM_obj_tol && change_proj < EAM_tol && feas_points>num_feas && abs(opt_dagger - q.'*theta_hash) > 1e-4 && abs(maxviol_hash) <EAM_maxviol_tol)
        theta_hat     = theta_hash;
        theta_optbound= opt_val;
        [c,CV] = KMS_31_Estep(theta_hash.',f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
        EI =  EI(1);
        flag_opt =1;
        return;
    end
    
    % Step 5) Update contraction counter
    if abs(change_proj) < 1e-6 || isnan(opt_EI_proj)
       % If change_proj = 0, contract parameter space
       contraction_counter = contraction_counter+1;
    elseif abs(opt_dagger - q.'*theta_hash) < 1e-4 && contraction_counter ~=0
       % If theta# is close to the boundary and we updated theta (change_proj > 1e-5), then expand the parameter space.
       contraction_counter = contraction_counter-1;
    end
    
    % Step 6) Update optimal value
    opt_val_old = opt_val;
end

% If failed to converge, output failure flag
theta_hat = theta_hash;
theta_optbound = opt_val;
[c,CV] = KMS_31_Estep(theta_hash.',f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions);
EI =  EI(1);
flag_opt = 0;

end








