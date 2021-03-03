function S = KMS_AUX2_drawpoints(e_point,dim_p,LB,UB,KMSoptions,A,b,x0)
%% Code Description: Draw points
% This function uses either Latin hypercube sampling or hit-and-run sampling
% to draw e_points from the parameter space [theta_lo,theta_hi] or
% {x : Ax <= b}.  
%
% INPUT:
%       e_point     Number of points to draw from hypercube
%       dim_p       Dimension of parameter space
%       LB          Lower bound on parameter space
%       UB          Upper bound on parameter space
%       A,b         Points satisfy Ax <= b.
%       
% OUTPUT:
%       S           e_point-by-dim_p matrix of points in the parameter space 

HR = KMSoptions.HR;
flag_error = 0;
if nargin <6
    %% Latin hypercube sampling
    S = zeros(e_point,dim_p);
    ran=rand(e_point,dim_p);

    % Construct S using uniform draws and permutations of
    % (0,1,2,...,e_point-1).
    for jj=1:dim_p
       S(:,jj) = (ran(:,jj) + (randperm(e_point) - 1).')/e_point;
    end

    S = repmat((UB-LB).',[e_point,1]).*S+repmat(LB.',[e_point,1]);
elseif HR == 0
    % If HR == 0, then we use uniform draw-and-discard sampling. 
    %This is a user-specified option.
    % This method draws points using latin hypercube sampling and checks
    % whether or not the satisfy the polytope constraints.
    % We keep drawing until we obtain e_points.  the
    flag_conv = 0;
    num_draw = e_point*100;
    S = [];
    while flag_conv == 0
        S_temp  = zeros(num_draw,dim_p);
        ran     = rand(num_draw,dim_p);

        % Construct S_temp using uniform draws and permutations of
        % (0,1,2,...,e_point-1).
        for jj=1:dim_p
            S_temp(:,jj) = (ran(:,jj) + (randperm(num_draw) - 1).')/num_draw;
        end

        % Map S_temp into [LB,UB]
        S_temp = repmat((UB-LB).',[num_draw,1]).*S_temp+repmat(LB.',[num_draw,1]);

        % Discard points not feasible
        ind = find(max(A*(S_temp.') - repmat(b,[1,size(S_temp,1)])) > 0).';
        S_temp(ind,:) = [];

        % Update S and terminate if we have a sufficient number of
        % points
        S = [S;S_temp];
        if size(S,1) >= e_point
            flag_conv = 1;
            S(e_point+1:end,:) = [];
        end
    end
else
    %% Hit-and-run sampling
    % If HR == 1, then we use hit-and-run sampling. This is a
    % user-specified option.
    A_aug = [A ; eye(size(UB,1)) ; -eye(size(LB,1))];
    b_aug = [b ; UB  ; -LB];
    if nargin <8
        options.method = 'achr';
        options.isotropic = 0;
        try
            S = cprnd(5*e_point,A_aug,b_aug,options);
        catch
            flag_error = 1;
        end
    else
       options.x0 = x0;
       options.method = 'achr';
       options.isotropic = 0;
       try
            S = cprnd(5*e_point,A_aug,b_aug,options);
       catch
            options.x0 = [];
            try 
                S = cprnd(5*e_point,A_aug,b_aug,options);
            catch
                flag_error = 1;
            end
        end
    end
    if ~flag_error 
        % Drop S that is not unique
        S =  uniquetol(S,1e-6,'ByRows',true);

        % Double check points satisfy constraints

        ind = find(max(A_aug*(S.') - repmat(b_aug,[1,size(S,1)])) > 0).';
        S(ind,:) = [];
        S(e_point+1:end,:) = [];

        % Check to see if HAR sampled correctly.  We should have e_point
        % uniformly sampled points
        if size(S,1) < e_point
            flag_error =1;
        end
    end
end
   
if flag_error
   % HAR did not uniformly sample correctly.  We revert to vertex sampling.
   
   % Find vertices
   try
        V = lcon2vert(A_aug,b_aug);
   catch
      S = [];
      return;
   end
   if isempty(V)
       S = [];
       return;
   end
   
   % Number of vertices
   nv = size(V,1);
   
   % Draw weights
   weight = rand(nv,e_point);
   weight = weight./repmat(sum(weight),[nv,1]);
   weight = weight.';
   
   % Construct randomly sampled points
   S = weight*V;
   
   % Double check points satisfy constraints
   ind = find(max(A_aug*(S.') - repmat(b_aug,[1,size(S,1)])) > 0).';
   S(ind,:) = [];
   S(e_point+1:end,:) = [];
end



end