%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Defines the GMS penalization function for DR inference as in Eq. (2.11)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = phi_function(xi,p,v)
% xi is the xi(1) vector in AS (2010) Eq. 4.5;
% p is the number of moment inequalities;
% v is the number of moment equalities;

% My definition of "infinity" (satisfies Inf*0 = 0 (as desired), while the built in "inf" has inf*0 = NaN)
Inf = 10^10;

% define GMS penalization;
value      = zeros(1,p+v); % zero for equalities and "close" inequalities;
value(1:p) = (xi(:,1:p)>1).*Inf; % infinity for "sufficiently violated" inequalities;
end