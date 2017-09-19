%This function implements a "hybrid" test for a non-linear parameter using linear programming

%The moments are assumed to be of the form y_T - X_T(lambda) * delta <= 0

%This function uses a hybrid approach to test whether lambda is in the
%identified set, where the LF approach is applied with size gamma, and then
%the conditional linear programming appraoch is applied with size
% 1 - (1-gamma)(1-alpha)

%Its arguments are:

%y_T : used in definition of moments; see above 
%X_T : used in definition of moments; see above
%Sigma: conditional variance of y_T | X_T. Assumed to be normalized to have
    %diagonal of ones
%gamma: signficance for LF test
% eta_vec (optional): the draws used when calculating the least favorable
% critical values. If not supplied, this is simulated inside the function using 1000 draws. 
    %If, however, the LF critical values are being calculated elsewhere, it may be
    %computationally more efficient to pass these draws here; users can also customize the 
    %number of draws by doing it outside

function reject = lp_hybrid_test_fn( y_T, X_T, Sigma, alpha, gamma, ds, varargin)

%% Compute the LF cutoff, if it's not given
%Set c_gamma if additional argument is provided; otherwise, compute it
if( isempty(varargin) == 0)
    eta_vec = varargin{1};
    c_gamma = quantile( eta_vec, 1-alpha);
else
    rng(0);
    numsims_lp = 1000;
    Z_draws = randn( size(y_T,1) , numsims_lp);
    c_gamma = c_lf_lp(X_T,Z_draws,Sigma,gamma);
end

%% Run the LP to calculate eta
%Store number of parameters and moments
M = size(Sigma,1);
k = size(X_T, 2);

%Compute eta, and the argmin delta
[eta, delta, lambda,error_flag] = test_delta_lp_fn( y_T, X_T, optimoptions('linprog','Algorithm','dual-simplex','TolFun', 10^-8, 'Display', 'off', 'MaxIter', 100000));

if(error_flag > 0)
    reject = 0;
    warning('LP for eta did not converge properly. Not rejecting');
    return;
end


%% Compare eta to the LF cutoff
%Reject if eta is greater than the LF cutoff and return
if(eta > c_gamma)
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
X_TB = X_T(B_index,:);
fullrank = rank(X_TB) == min( size(X_TB) );



if(~fullrank || degenerate)
    
    %Calculate vlo and vup using the bisection approach that conditions on
    %having a gamma_tilde - a vertex of the dual (note that since matlab
    %implements the dual-simplex method, lambda is guaranteed to be such a
    %gamme_tilde
    [vlo_dual,vup_dual,eta_dual,gamma_tilde] = lp_dual_fn( y_T, X_T, eta,lambda, Sigma);

    vup_dual = min([vup_dual, c_gamma]); %this conditions on having not rejected the LF test 
    
    sigma_B_dual = sqrt( gamma_tilde' * Sigma * gamma_tilde);
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
        alpha_tilde = (alpha - gamma) / (1 - gamma);
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

X_TB = X_T(B_index,:);
X_TBc = X_T(Bc_index,:);

S_B = eye(M);
S_B = S_B(B_index, :);

S_Bc = eye(M);
S_Bc = S_Bc(Bc_index, :);


Gamma_B  = [i_Bc , X_TBc] * [i_B , X_TB]^(-1) * S_B - S_Bc;

e1 = [1 ; zeros(size_B - 1 ,1) ];
v_B = ( e1' * [i_B , X_TB]^(-1) * S_B )';

Gamma_B_tilde = [Gamma_B; -v_B'];
u = [ zeros( size( Gamma_B, 1), 1) ; -c_gamma];

sigma2_B = v_B' * Sigma * v_B;
sigma_B = sqrt(sigma2_B);

rho = Gamma_B_tilde * Sigma * v_B / sigma2_B;


maximand_or_minimand = (u -Gamma_B_tilde * y_T )./ rho + v_B' * y_T;

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
    alpha_tilde = (alpha - gamma) / (1 - gamma);
    reject = pval < alpha_tilde;

end



end
