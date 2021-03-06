function [reject, eta] = conditional_test_fn( y, X, Sigma, alpha)

%% Compute the conditional test from Andrews, Roth and Pakes for testing moments of the form E[Y_i - X_i delta | Z_i] <= 0
%Inputs
% y: the (scaled) sample average of Y_i (a k x 1 vector)
% X: the (scaled) sample average of X_i (a k x m vector)
% Sigma: the (scaled) estimate of E[ Var(Y_i|Z_i) ] 
% alpha: the size of the test (e.g. 0.05 for 5% significance)
%Outputs:
% reject: does the test reject
% eta: the test statistic
%Notes
% If y and X are sample averages, Sigma should be E[ Var(Y|X) ] /N
% If Sigma is E[ Var(Y|X) ], then y and X should be averages scaled by sqrt(N)


%Check if Sigma has 1s on the diagnol. If not, renormalize. (Original code
%assumes Sigma is a correlation matrix)
if( max(abs(diag(Sigma)-1)) > 10^-6 )
    Sigma_invsqrt = diag( diag( Sigma^(-1/2) ) );
    [reject, eta] = conditional_test_fn(Sigma_invsqrt * y, ...
                                           Sigma_invsqrt * X, ...
                                           Sigma_invsqrt * Sigma * Sigma_invsqrt',...
                                           alpha );
    return;
end
%Store number of parameters and moments
M = size(Sigma,1);
k = size(X, 2);
%Compute eta, and the argmin delta

%[eta, delta, lambda] = etahat_fn( y_T, X_T, optimoptions('linprog','Algorithm','interior-point'));


%tol = 10^(-6);
%B_index = abs(lambda) > tol; %the bidning moments are the ones at which the Lagrange multiplier is >tol
%Bc_index = B_index == 0;


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


%%The following block checks for conditions under which the primal
%%%and dual solutions are equal to one another. If these conditions don't hold, we go to the dual
 
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


if( (~fullrank) || degenerate)
    warning('Using bisection approach for vlo/vup');
    %Calculate vlo and vup using the bisection approach that conditions on
    %having a gamma_tilde - a vertex of the dual (note that since matlab
    %implements the dual-simplex method, lambda is guaranteed to be such a
    %gamme_tilde
    [vlo_dual,vup_dual,eta_dual,gamma_tilde] = lp_dual_fn( y, X, eta,lambda, Sigma);
    
    sigma_B_dual = sqrt( gamma_tilde' * Sigma * gamma_tilde);
    maxstat = eta_dual ./ sigma_B_dual;
    zlo_dual = vlo_dual ./ sigma_B_dual;
    zup_dual = vup_dual ./ sigma_B_dual;
    
    
    if( ~ (zlo_dual <= maxstat && maxstat <= zup_dual) )
        warning(strcat('max_stat (', num2str(maxstat), ')',...
                   'is not between z_lo (', num2str(zlo_dual), ')',...
                   'and z_up (', num2str(zup_dual), ...
                   ') in the dual approach in dataset ', num2str(ds)) );
    
        reject = 0;
        return;
    else
        pval = Truncated_normal_p_value_by_simulation(maxstat,zlo_dual,zup_dual);
        reject = pval < alpha;
    return;
    end


 
end

%warning('Things look good. Using primal');

 %
 
%  slack = y_T - X_T * delta - eta;
%  [~, sorted_index] = sort(slack, 'descend');
%  smallest_kplus1_moments = sorted_index(1:(k+1));
%  B_index = false(size(slack));
%  B_index(smallest_kplus1_moments) = true;
%  Bc_index = B_index == 0;

%Compute ingredients for the test
size_B = sum(B_index == 1);
size_Bc = sum(B_index == 0);

i_B = ones(size_B,1);
i_Bc = ones(size_Bc,1);

%X_TB = X_T(B_index,:); %this is now done earlier
X_TBc = X(Bc_index,:);

S_B = eye(M);
S_B = S_B(B_index, :);

S_Bc = eye(M);
S_Bc = S_Bc(Bc_index, :);


Gamma_B  = [i_Bc , X_TBc] * [i_B , X_TB]^(-1) * S_B - S_Bc;

e1 = [1 ; zeros(size_B - 1 ,1) ];

v_B = ( e1' * [i_B , X_TB]^(-1) * S_B )';
sigma2_B = v_B' * Sigma * v_B;
sigma_B = sqrt(sigma2_B);

rho = Gamma_B * Sigma * v_B / sigma2_B;


maximand_or_minimand = -Gamma_B * y ./ rho + v_B' * y;

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


z_lo = (v_lo / sigma_B);
z_up = (v_up / sigma_B);
max_stat = (eta / sigma_B);

if( ~ (z_lo <= max_stat && max_stat <= z_up) )
    warning(strcat('max_stat (', num2str(max_stat), ')',...
                   'is not between z_lo (', num2str(z_lo), ')',...
                   'and z_up (', num2str(z_up), ...
                   ') in the primal approach ') );
    reject = 0;

else
    pval = Truncated_normal_p_value_by_simulation( max_stat, z_lo , z_up );
    reject = pval < alpha;

end



end
