%This function computes vlo and vup for the dual linear program for the
%conditional values

%Inputs: see "Conditioning via the Dual" note for notation
function [vlo, vup] = vlo_vup_dual_fn(eta,s_T,gamma_tilde, Sigma, W_T)

tol_c = 10^(-6);
tol_equality = 10^(-6);

sigma_B = sqrt( gamma_tilde' * Sigma * gamma_tilde);

lowinitial = min(-100,eta - 20 * sigma_B);
highinitial = max(100,eta + 20 * sigma_B);
maxiters = 10000;


%% Compute vup
dif = tol_c +1;
iters = 0;


low = eta;
high = highinitial;


check_if_solution = @(c) check_if_solution_helper(c, tol_equality, s_T, gamma_tilde,Sigma,W_T);

if(~ check_if_solution(eta) )
    warning('The user-supplied eta does not appear to be a solution. Not rejecting automatically');
    vlo = eta; 
    vup = Inf;
    return;
end

if( check_if_solution( highinitial) )
    vup = inf;
else
    
    %Throughout the while loop, the high value is not a solution and the
    %low-value is a solution
    
    %If the midpoint between them is a solution, then we set low to the
    %midpoint; otherwise, we set
    while( dif>tol_c && iters < maxiters)
        
        c = 1/2 * (low + high);
        
        
        if( check_if_solution(c) )
            low = c;
        else
            high = c;
        end
        
        
        
        dif = high - low;
        iters = iters + 1;
        
    end
    
    
    if(iters == maxiters)
        warnings('Reached max iters without a solution');
    end
    
    vup = c;
end



%% Compute vlo
dif = tol_c +1;
iters = 0;


low = lowinitial;
high = eta;


if( check_if_solution( lowinitial) )
    vlo = -inf;
else
    
    %Throughout the while loop, the low value is not a solution and the
    %high-value is a solution
    
    %If the midpoint between them is a solution, then we set high to the
    %midpoint; otherwise, we set low to the midpoint
    while( dif>tol_c && iters < maxiters)
        
        c = 1/2 * (low + high);
        
        
        if( check_if_solution(c) )
            high = c;
        else
            low = c;
        end
        
        
        
        dif = high - low;
        iters = iters + 1;
        
    end
    
    
    if(iters == maxiters)
        warnings('Reached max iters without a solution');
    end
    
    vlo = c;
end



end


function maxvalue = max_program( s_T, gamma_tilde,Sigma,W_T, c)

% c = min_gamma (s_T + Sigma*gamma_tilde / (gamma_tilde' * Sigma *
% gamma_tilde) * c)' * gamma

% s.t

% W_T' * gamma = e_1
% gamma >= 0

f = (s_T + Sigma * gamma_tilde ./ (gamma_tilde' * Sigma * gamma_tilde) .* c);
Aeq = W_T';
beq = [1; zeros( size(Aeq,1) - 1, 1) ]; 
lb = zeros( size(f) );

[~,maxvalue] = linprog(-f,[],[],Aeq,beq,lb,[],[], optimoptions('linprog','TolFun', 10^(-8), 'Display','off') );
maxvalue = -maxvalue; %Because we minimize the negative to get the max
end

function equal = check_equals(a ,b , tol)
        equal = abs(a-b) <= tol;
end

function solution = check_if_solution_helper(c, tol, s_T, gamma_tilde,Sigma,W_T)
    
    minvalue = max_program( s_T, gamma_tilde,Sigma,W_T, c);
    solution = check_equals(c,minvalue,tol);
end