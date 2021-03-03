function c_val = KMS_32_Critval(phi_test,f_ineq_keep,f_eq_keep,G_ineq,G_eq,Dg_ineq,Dg_eq,A_rho,b_rho,KMSoptions)
%% Critical value
% This function computes the critical value for theta_test.
% The critical value is the fixed point to the following problem:
%
%       h(c) = (1/B) sum_b psi_b(c) - (1-alpha)
%
% That is, c_val is the real number satisfying h(c_val) = 0.
% We compute c_val using Brent's fixed-point method.  Note that, we take
% advantage of the property that psi_b is monotone in c and
% that psi_b is a binary variable.  That is, we do not need to recompute
% psi_b(c') for psi_b(c) = 1 if c' >c.  Hence the reason why we do not use the
% optimized fzero MatLab command.
%
% Also, note that (1/B) sum_b psi_b(.) takes on values
% 0,(1/B),(2/B),...,(B-1)/B,1.  So the best that we can do is get within
% (1/2B) of 0.  So, we set the tolerance equal to 1/B.
%
% For more information on Brent's fixed-point algorithm see
%       https://en.wikipedia.org/wiki/Brent%27s_method
% or Numerical Analysis by Burden, Faires, Burden, or any other good source
% on numerical analysis.
%
% INPUTS:
%   phi_test            GMS function evaluated at theta_test
%
%   f_ineq_keep,f_eq_keep   Moments to keep   
%  
%   G_ineq,G_eq         Bootstrapped and recentered moments
%
%   Dg_ineq,Dg_eq       Gradients of moment function evaluated at theta_test
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
%   c_val               If KMS confidence interval requested, c_val is the 
%                       critical value evaluated at theta_test, i.e., the
%                       real number that solves: 
%                       h(c_val) = (1/B) sum_b psi_b(c) - (1-alpha) = 0
%                       If AS confidence interval is requested, c_val is
%                       the (1-alpha) quantile of the Andrew and Soares
%                       bootstrapped distribution: max(G+phi_test)_+

%% Extract relevant information from KMSoptions
c_B         = KMSoptions.c_B;
B           = KMSoptions.B;
alpha       = KMSoptions.alpha;
c_maxit     = KMSoptions.c_maxit;                          % Maximum number of iterations in root-finding algorithm for critical val
hc_tol      = KMSoptions.hc_tol;                           % Tolerence level for root-finding algorithm
c_tol       = KMSoptions.c_tol;                            % Tolerence level for root-finding algorithm
Brent_tol   = KMSoptions.Brent_tol;                        % Tolerence level for root-finding algorithm -- if difference is close to
                                                           % zero, then use Bisection method.
CI_method = KMSoptions.CI_method;                          % Confidence interval method

%% Calculate the critical value
if strcmp(CI_method,'KMS')
%% Intialize
% Set interval to be [a,b] = [0,c_B], where c_B is the Bonferroni critical value.
% We will find c on [a,b] such that h(c) = 0.
a       = 0;
b       = c_B;
psi_a   = KMS_33_Coverage(a,zeros(B,1),ones(B,1),phi_test,f_ineq_keep,f_eq_keep,G_ineq,G_eq,Dg_ineq,Dg_eq,A_rho,b_rho,KMSoptions);
h_a     = mean(psi_a) - (1-alpha);
if h_a > 0 
   % If theta is deep inside the constraint set it is possible that the
   % critical value is c_val = 0.
   c_val = 0;
   return
end
psi_b   = KMS_33_Coverage(b,psi_a,ones(B,1),phi_test,f_ineq_keep,f_eq_keep,G_ineq,G_eq,Dg_ineq,Dg_eq,A_rho,b_rho,KMSoptions);
h_b     = mean(psi_b) - (1-alpha);

% Need to make sure h_b > 0.  If not, make b larger.
if h_b < 0
   for jj = 1:1000
       b = 2*b;
       psi_b   = KMS_33_Coverage(b,psi_a,ones(B,1),phi_test,f_ineq_keep,f_eq_keep,G_ineq,G_eq,Dg_ineq,Dg_eq,A_rho,b_rho,KMSoptions);
       h_b     = mean(psi_b) - (1-alpha);
       if h_b > 0
           break;
       end
   end
   if jj == 1000 
      error('Cannot find root.  Possibly a problem with gradient function. Check to make sure gradients are not NaN or inf here.')
   end
end

% Note: psi(c) is a vector of indicators of whether or not the bth
% bootstrap repetition evaluated at (c) intersects the hyperplane
% p'lambda=0.

% Note: The vector psi_n is updated based c_{n} and c_{n-1} for
% computational reasons:
%    Case 1:  If c_b > c_a, then psi_b >= psi_a.
%             So if psi_a=1, then we do not need to recompute
%             psi_b -- it is equal to 1.  Therefore we only
%             update psi_b such that psi_a = 0
%   Case 2:   If c_b < c_a, then psi_b <= psi_a.
%             So if psi_a=0, then we do not need to recompute
%             psi_b -- it is equal to 0.  Therefore we only
%             update psi_b such that psi_a = 1

% Run root finding algorithm only if h_b not equal to zero (it could be the
% case that c_B is the root).
if abs(h_b) < c_tol
    c_val = b;
    return;
end

% Initiate by setting c = a
c = a;
h_c = h_a;

% Set flag
% This determines if we use bisection method, quadratic, or secant method.
mflag = 1;
for ii = 1:c_maxit
    if h_a ~= h_c && h_b ~= h_c
        % Inverse quadratic interpolation
        s = a*h_b*h_c/((h_a - h_b)*(h_a - h_c)) ...
            + b*h_a*h_c/((h_b  - h_a)*(h_b - h_c)) ...
            + c*h_a*h_b/((h_c - h_a)*(h_c - h_b));
    else
        % Secant method
        s = b - h_b*(b-a)/(h_b -h_a);
    end

    % Check conditions to determine if we should use s above or
    % bisection method
    if (s < (3*a+b)/4 || s > b) ...
            || (mflag == 1 && 2*abs(s-b) >= abs(b-c)) ...
            || (mflag == 0 && 2*abs(s-b) >= abs(c-d)) ...
            || (mflag == 1 && Brent_tol > abs(b-c)) ...
            || (mflag == 0 && Brent_tol > abs(c-d))
        % Use bisection, update flag
        s = 0.5*(a+b);
        mflag = 1;
    else
        % Use sectant/quadratic, set flag to 0.
        mflag = 0;
    end

    % Calculate h(s).  To do this, first update psi(s) using the
    % best guess so far, which is (b,psi(b)):
    psi_s = KMS_33_Coverage(s,psi_a,psi_b,phi_test,f_ineq_keep,f_eq_keep,G_ineq,G_eq,Dg_ineq,Dg_eq,A_rho,b_rho,KMSoptions);
    h_s = mean(psi_s) - (1-alpha);

    % Update d and c.  Note that d tracks the value of c from iteration k-2
    % d is used to check which method we should use.
    d = c;
    c = b;
    h_c = h_b;

    if h_a*h_s < 0
        % Root is bracketed by [a,s], so update by setting b=s
        b = s;
        psi_b = psi_s;
        h_b = h_s;
    else
        % Otherwise root is bracketed by [s,b], so update by setting a=s
        a = s;
        psi_a = psi_s;
        h_a = h_s;
    end

    % Check ordering of h_a,h_b.  We always set b to be our "best root".
    if abs(h_a) < abs(h_b)
        % Swap (a,b), since h_a is better guess in this case.
        x = a;
        psi_x = psi_a;
        h_x = h_a;
        a = b;
        psi_a = psi_b;
        h_a = h_b;
        b = x;
        psi_b = psi_x;
        h_b = h_x;
    end

    % Tolerence check
    % We only output b iff we are close to optimal AND coverage >=
    % (1-alpha).  We do not output if coverage < 1-alpha, because in the
    % worst case scenario we would rather have conservative bias.
    if (abs(h_b) < hc_tol || abs(b-a) < c_tol) && h_b >= 0
        c_val = b;
        return;
    end
end
% If failed to converge in the number of iterations, output b if h_b >=0,
% otherwise output a.  In this case, the KMS interval will suffer from
% conservative bias.
if h_b >= 0
    c_val = b;
else
    c_val = a;
end

elseif strcmp(CI_method,'AS')
G = [G_ineq;G_eq];
B = size(G,2);
S_3_boot = max(0,max(G+repmat(phi_test,1,B)));
c_val = quantile(S_3_boot,1-alpha);
end


end
