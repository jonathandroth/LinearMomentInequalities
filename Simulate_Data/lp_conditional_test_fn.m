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
[eta, delta, lambda,flag] = test_delta_lp_fn( y_T, X_T, optimoptions('linprog','Algorithm','dual-simplex', 'Display', 'off'));

if(flag >0 )
    error('Trying to do conditional test with infinite cutoff');

end

%%Store which moments are binding
    %Note we manually calculate this (within a tolerance), rather than
    %using the lagrange multipliers, since these can sometimes be 0 at a
    %binding moments
 tol = 10^(-6);
 slack = y_T - X_T * delta - eta;
 B_index = abs(slack) < tol;
 Bc_index = B_index == 0;
 
 
 %Check if the right number of moments are binding. If not, throw a warning
 %and return NA
 if( sum(B_index) ~= (k+1) )
     if( sum(B_index) > (k+1) )
        warning('Number of Binding Moments More Than k+1');
        %reject = NaN;
        reject = 0;
        return;
     else
        warning('Number of Binding Moments Less Than k+1');
        %reject = NaN;
        reject = 0;
        
        return;
     end
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

X_TB = X_T(B_index,:);
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
zeta_lo = normcdf( v_lo / sigma_B);
zeta_up = normcdf( v_up / sigma_B);


reject = (eta / sigma_B) >= norminv( (1- alpha) * zeta_up + alpha * zeta_lo );

if( norminv( (1- alpha) * zeta_up + alpha * zeta_lo ) == Inf)
    warning('Infinite Critical Value Computed');
    
end
end