
function [reject, eta] = hybrid_test_fn( y, X, Sigma, alpha, kappa, varargin)
%% Computing the hybrid test from Andrews, Roth and Pakes for testing moments of the form E[Y_i - X_i delta | Z_i] <= 0
%Inputs:
% y: the (scaled) sample average of Y_i (a k x 1 vector)
% X: the (scaled) sample average of X_i (a k x m vector)
% Sigma: the (scaled) estimate of E[ Var(Y_i|Z_i) ] 
% alpha: the size of the test (e.g. 0.05 for 5% significance)
% kappa: the size of the first-stage LF test; we recommend alpha/10
% eta_vec (optional): provide the eta_vec output from the
% lf_critical_value_fn to increase speed of the hybrid test
%Outputs:
% reject: does the test reject
% eta: the test statistic
%Notes:
% If y and X are sample averages, Sigma should be E[ Var(Y|X) ] /N
% If Sigma is E[ Var(Y|X) ], then y and X should be averages scaled by sqrt(N)



%Check if Sigma has 1s on the diagnol. If not, renormalize. (Original code
%assumes Sigma is a correlation matrix)
if( max(abs(diag(Sigma)-1)) > 10^-6 )
    Sigma_invsqrt = diag( diag( Sigma^(-1/2) ) );
    
    if(isempty(varargin))
        [reject, eta] = hybrid_test_fn(Sigma_invsqrt * y, ...
                                               Sigma_invsqrt * X, ...
                                               Sigma_invsqrt * Sigma * Sigma_invsqrt',...
                                               alpha, kappa);
        return;

    else
        [reject, eta] = hybrid_test_fn(Sigma_invsqrt * y, ...
            Sigma_invsqrt * X, ...
            Sigma_invsqrt * Sigma * Sigma_invsqrt',...
            alpha, kappa, varargin{1});
        return;
    end     
        
end

%% Run the LP to calculate eta
%Store number of parameters and moments
M = size(Sigma,1);
k = size(X, 2);

%Compute eta, and the argmin delta
[eta, delta, lambda,error_flag] = etahat_fn( y, X, Sigma, optimoptions('linprog','Algorithm','dual-simplex','TolFun', 10^-8, 'Display', 'off', 'MaxIter', 100000));

if(error_flag > 0)
    reject = 0;
    warning('LP for eta did not converge properly. Not rejecting');
    return;
end


%Check if eta<0. If so, don't reject
if(eta <0)
    reject =0;
    return;
end    

%% Compute the LF cutoff, if it's not given
%Set c_kappa if additional argument is provided; otherwise, compute it
if( isempty(varargin) == 0)
    eta_vec = varargin{1};
    c_kappa = quantile( eta_vec, 1-alpha);
else
    rng(0);
    numsims_lp = 1000;
    Z_draws = randn( size(y,1) , numsims_lp);
    c_kappa = lf_critical_value_fn(X,Z_draws,Sigma,kappa);
end



%% Compare eta to the LF cutoff
%Reject if eta is greater than the LF cutoff and return
if(eta > c_kappa)
    reject = 1;
    return;
end

%% If don't reject with LF, do the conditional test

 
%Check whether problem is degenerate 
tol_lambda = 10^(-6);
degenerate = sum( lambda>tol_lambda ) ~= (k+1) ;


%%Store which moments are binding
%%%Currently doing this using lambda rather than moments to avoid differing
%%%precision issues
 %tol_slack = 10^(-6);
 %slack = y_T - X_T * delta - eta;
 %B_index = abs(slack) < tol_slack;
  
 B_index = lambda > tol_lambda;
 Bc_index = B_index == 0;
 
%Check whether X_T,B has full rank
X_TB = X(B_index,:);
fullrank = rank(X_TB) == min( size(X_TB) );



if(~fullrank || degenerate)
     warning('Using bisection approach for vlo/vup');
    %Calculate vlo and vup using the bisection approach that conditions on
    %having a kappa_tilde - a vertex of the dual (note that since matlab
    %implements the dual-simplex method, lambda is guaranteed to be such a
    %gamme_tilde
    [vlo_dual,vup_dual,eta_dual,kappa_tilde] = lp_dual_fn( y, X, eta,lambda, Sigma);

    vup_dual = min([vup_dual, c_kappa]); %this conditions on having not rejected the LF test 
    
    sigma_B_dual = sqrt( kappa_tilde' * Sigma * kappa_tilde);
    max_stat = eta_dual ./ sigma_B_dual;
    zlo_dual = vlo_dual ./ sigma_B_dual;
    zup_dual = vup_dual ./ sigma_B_dual;
    
    
    if( ~ (zlo_dual <= max_stat && max_stat <= zup_dual) )
        warning(strcat('Hybrid: max_stat (', num2str(max_stat), ')',...
                   'is not between z_lo (', num2str(zlo_dual), ')',...
                   'and z_up (', num2str(zup_dual), ...
                   ') in the dual approach in dataset ', num2str(ds)) );
    
        reject = 0;
        return;
    else
        pval = Truncated_normal_p_value_by_simulation(max_stat,zlo_dual,zup_dual);
        alpha_tilde = (alpha - kappa) / (1 - kappa);
        reject = pval < alpha_tilde;
    return;
    end
    
    
    
end

%warning('Things look good. Using primal')

%% If not degenerate, then use the "primal approach"
%Compute ingredients for the test
size_B = sum(B_index == 1);
size_Bc = sum(B_index == 0);

i_B = ones(size_B,1);
i_Bc = ones(size_Bc,1);

X_TB = X(B_index,:);
X_TBc = X(Bc_index,:);

S_B = eye(M);
S_B = S_B(B_index, :);

S_Bc = eye(M);
S_Bc = S_Bc(Bc_index, :);


kappa_B  = [i_Bc , X_TBc] * [i_B , X_TB]^(-1) * S_B - S_Bc;

e1 = [1 ; zeros(size_B - 1 ,1) ];
v_B = ( e1' * [i_B , X_TB]^(-1) * S_B )';

kappa_B_tilde = [kappa_B; -v_B'];
u = [ zeros( size( kappa_B, 1), 1) ; -c_kappa];

sigma2_B = v_B' * Sigma * v_B;
sigma_B = sqrt(sigma2_B);

rho = kappa_B_tilde * Sigma * v_B / sigma2_B;


maximand_or_minimand = (u -kappa_B_tilde * y )./ rho + v_B' * y;

if( sum(rho >0) > 0 )
    v_lo = max( maximand_or_minimand( rho >0) );
else
    v_lo = -Inf;
end

if( sum(rho <0) >0 )
    v_up = min( maximand_or_minimand( rho <0) );
else
    v_up = Inf;
end

z_lo = v_lo / sigma_B; %do i need to add c_0 here?
z_up = v_up / sigma_B;%do i need to add c_0 here?
max_stat = eta ./ sigma_B;


if( ~ (z_lo <= max_stat && max_stat <= z_up) )
    warning(strcat('Hybrid: max_stat (', num2str(max_stat), ')',...
                   'is not between z_lo (', num2str(z_lo), ')',...
                   'and z_up (', num2str(z_up), ...
                   ') in the primal approach in dataset ', num2str(ds)) );
    reject = 0;

else
    pval = Truncated_normal_p_value_by_simulation(max_stat,z_lo,z_up);
    alpha_tilde = (alpha - kappa) / (1 - kappa);
    reject = pval < alpha_tilde;

end



end
