function [theta_keep,EI_keep]  = KMS_36_drawpoints(theta_hash,q,r_max,r_min,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,dmodel,LB,UB,A,b,KMSoptions)
%% Code description
% This code draws initial points theta near theta_hash.  It checks that at
% these draws satisfy EI>0.
%
% Method 1: Draw points from a box:
% A  = { theta : p'theta >= theta_hash, d(theta,theta_hash) < r}
% for some tuning parameter r, where the distance metric in the infinity
% norm.  We start with r = r_max and add, say, 2 points and check expected
% improvement of these two points.  If EI > 0 we add them, otherwise
% we continue.  We reduce the size of the box: r = r/2, and add more
% points. We continue until we have either, say, 10 points or r<r_min.
%
% Method 2: Draw points from a slice:
% Draw points from
% B  = { theta : p'theta >= theta_hash, and |theta(component) -theta_hash(component)| < r}.
% Again, we continue drawing until we have a
% sufficient number of points satisfying EI>0.
%
% Method 3: Draw some points uniformly from the parameter space

%% Extract Parameters
dim_p           = KMSoptions.dim_p;
component       = KMSoptions.component;
unif_points     = KMSoptions.unif_points;
EI_points       = KMSoptions.EI_points;
EI_points_start = KMSoptions.EI_points_start;
parallel        = KMSoptions.parallel;
sample_method   = KMSoptions.sample_method;
options_linprog = KMSoptions.options_linprog;

% Check to make sure r_max is biggest
r_max = max(r_max, r_min);

% Initiate:
r = r_max;
theta_keep = [];
EI_keep    = [];

LB2 = LB;
UB2 = UB;
A1 = A;
b1 = b;
A2 = A;
b2 = b;

%% Draw points
flag_while = true;
while flag_while
    if sample_method == 0   % Latin hypercube sampling
        % Method 1: Draw points from a box around theta#
         if max(q) > 0
            LB1 = LB;
            UB1 = min(theta_hash + r,UB);
        else
            LB1 = max(theta_hash - r,LB);
            UB1 = UB;
        end
        theta_draw1 = KMS_AUX2_drawpoints(EI_points_start,dim_p,LB1,UB1,KMSoptions);
        
        % Method 2: Draw points from slice around theta#
        if max(q) > 0
            LB2 = LB;
            UB2 = UB;
            UB2(component) = min(theta_hash(component) + r,UB(component));
        else
            LB2 = LB;
            LB2(component) = max(theta_hash(component) - r,LB(component));
            UB2 = UB;
        end
        theta_draw2 = KMS_AUX2_drawpoints(EI_points_start,dim_p,LB2,UB2,KMSoptions);
    
    elseif sample_method == 1 || sample_method == 2
        % Hit-and-run sampling
        % Method 1: Draw points from a box around theta#
        A1 = A;
        A1 = [A1; eye(dim_p); -eye(dim_p)];
        b1 = b;
        b1 = [b1 ; theta_hash+r; -theta_hash+r];
        theta_draw1 = KMS_AUX2_drawpoints(EI_points_start,dim_p,LB,UB,KMSoptions,A1,b1,theta_hash);
        
        % Method 2: Draw points from a slice around theta#
        A2 = A;
        A2 = [A2; q.'];
        b2 = b;
        b2 = [b2 ; q.'*theta_hash + r*abs(q.'*ones(dim_p,1))];
        theta_draw2 = KMS_AUX2_drawpoints(EI_points_start,dim_p,LB,UB,KMSoptions,A2,b2,theta_hash); 
    elseif sample_method == 3
        % Method 1: Draw points from a box around theta#
        if max(q) > 0
            LB1 = LB;
            UB1 = min(theta_hash + r,UB);
        else
            LB1 = max(theta_hash - r,LB);
            UB1 = UB;
        end
    
        % Bound transform if draw-and-discard method is used
        [LB1,UB1] = bound_transform(LB1,UB1,KMSoptions);
        
        theta_draw1 = KMS_AUX2_drawpoints(EI_points_start,dim_p,LB1,UB1,KMSoptions,A,b,theta_hash);
        
        % Method 2: Draw points from slice around theta#
        component = find(abs(q) ==1);
        if max(q) > 0
            LB2 = LB;
            UB2 = UB;
            UB2(component) = min(theta_hash(component) + r,UB(component));
        else
            LB2 = LB;
            LB2(component) = max(theta_hash(component) - r,LB(component));
            UB2 = UB;
        end
        
        % Bound transform if draw-and-discard method is used
        [LB2,UB2] = bound_transform(LB2,UB2,KMSoptions);
        
        theta_draw2 = KMS_AUX2_drawpoints(EI_points_start,dim_p,LB2,UB2,KMSoptions,A,b,theta_hash);
    end
    
    % Concatenate theta_draw1,draw2
    theta_draw = [theta_draw1;theta_draw2];
    size_draw = size(theta_draw,1);
    
    % Find theta's that have positive EI
    Eimprovement = @(theta)KMS_37_EI_value(theta,q,theta_hash,f_ineq,f_eq,f_ineq_keep,f_eq_keep,f_stdev_ineq,f_stdev_eq,dmodel,KMSoptions);
    EI = zeros(size_draw,1);
    if parallel
        parfor jj = 1:size_draw
            EI(jj,1) = -(max(Eimprovement( theta_draw(jj,:).')));
        end
    else
        for jj = 1:size_draw
            EI(jj,1) = -(max(Eimprovement( theta_draw(jj,:).')));
        end
    end
    
    % Sort theta from best, and pick those that have positive EI.
    [EI,I] = sort(EI,'descend');
    theta_draw = theta_draw(I,:);
    ind = find(EI>1e-10);
    
    % Keep those with positive EI
    theta_keep = [theta_keep ; theta_draw(ind,:)];
    EI_keep    = [EI_keep  ;EI(ind)];
    
    % Update r
    r = r/2;
    
    % If r < r_min or if we have sufficient number of points, break.
    if r < r_min || size(theta_keep,1)  > EI_points
        flag_while = false;
    end
end

% We also draw some points uniformly from
% {theta : p'theta >= p'theta#}
if sample_method == 0   % Latin hypercube sampling
    theta_draw = KMS_AUX2_drawpoints(unif_points,dim_p,LB,UB,KMSoptions);
else
    theta_draw = KMS_AUX2_drawpoints(unif_points,dim_p,LB,UB,KMSoptions,A,b,theta_hash);
end
% Set uniformly drawn EI to zero (this might not be true, but we do not
% require EI to be positive for these points so there is no need to
% calculate it.
EI = zeros(size(theta_draw,1),1);

% Concatenate
theta_keep = [theta_keep;theta_draw];
EI_keep  =[EI_keep ; EI];

end