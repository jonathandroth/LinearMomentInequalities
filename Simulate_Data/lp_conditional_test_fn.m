function reject = lp_conditional_test_fn( y_T, X_T, Sigma, alpha)

%Store number of parameters and moments
M = size(Sigma,1);
k = size(X_T, 2);
%Compute eta, and the argmin delta

%[eta, delta, lambda] = test_delta_lp_fn( y_T, X_T, optimoptions('linprog','Algorithm','interior-point'));


%tol = 10^(-6);
%B_index = abs(lambda) > tol; %the bidning moments are the ones at which the Lagrange multiplier is >tol
%Bc_index = B_index == 0;


%Compute eta, and the argmin delta
[eta, delta, lambda,error_flag] = test_delta_lp_fn( y_T, X_T, optimoptions('linprog','Algorithm','interior-point', 'Display', 'off', 'MaxIter', 100000));


if(error_flag >0 )
    error('Trying to do conditional test with infinite cutoff');

end

%%Store which moments are binding
 tol_slack = 10^(-6);
 slack = y_T - X_T * delta - eta;
 B_index = abs(slack) < tol_slack;
 Bc_index = B_index == 0;
 
%Check whether X_T,B has full rank
X_TB = X_T(B_index,:);
fullrank = rank(X_TB) == min( size(X_TB) );
 
%Check whether problem is degenerate 
tol_lambda = 10^(-6);
degenerate = sum( lambda>tol_lambda ) < (k+1) ;

if(~fullrank || degenerate)
    warning('Primal LP non-unique or degenerate. Using dual approach');
    [vlo_dual,vup_dual,eta_dual,gamma_tilde] = lp_dual_fn( y_T, X_T, Sigma);
    
    sigma_B_dual = sqrt( gamma_tilde' * Sigma * gamma_tilde);
    maxstat = eta_dual ./ sigma_B_dual;
    pval = Truncated_normal_p_value(maxstat,vlo_dual,vup_dual);
    
    reject = pval < alpha;
    return;
end

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
X_TBc = X_T(Bc_index,:);

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


maximand_or_minimand = -Gamma_B * y_T ./ rho + v_B' * y_T;

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

pval = Truncated_normal_p_value( (eta ./ sigma_B), v_lo, v_up);

reject = pval < alpha;

end
