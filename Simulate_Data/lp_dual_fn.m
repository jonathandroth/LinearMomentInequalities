function [vlo,vup,eta, gamma_tilde, error_in_lp] = lp_dual_fn( y_T, X_T, Sigma)


%% Solve the minmization for eta
W_T = [ ones(size(X_T,1) , 1) , X_T];

f = y_T;
Aeq = W_T';
beq = [1; zeros(size(Aeq,1)-1,1)];
lb = zeros( size(f) );

[gamma_tilde,eta, flag] = linprog(-f,[],[],Aeq,beq,lb,[],[], optimoptions('linprog','TolFun', 10^(-8), 'Display','off') );
eta = -eta; %Because we minimize the negative to get the max

if( flag ~= 1)
    warning('Error in dual approach linear program for eta. To be conservative, I do not reject');
    error_in_lp = 1;
else
    error_in_lp = 0;
end
%% Calculate vlo and vup
s_T =  ( eye( size(y_T,1) ) - (Sigma * (gamma_tilde * gamma_tilde')) ./ (gamma_tilde' * Sigma * gamma_tilde) ) * y_T;
[vlo, vup] = vlo_vup_dual_fn(eta,s_T,gamma_tilde, Sigma, W_T);

end