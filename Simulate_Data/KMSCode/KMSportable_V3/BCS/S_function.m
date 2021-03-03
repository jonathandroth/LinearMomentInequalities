%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Defines the S function according to MMM in Eq. (2.6);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Qn,chi] = S_function(m_std,p)
% m_std denotes the studentized sample moment conditions;
% p denotes the number of moment inequalities (which should appear first);

chi = 2; % degree of homogeneity of criterion function in MMM

% take negative part of inequalities;
m_std(:,1:p) = min(m_std(:,1:p),zeros(size(m_std,1),p));

% sample criterion function;
Qn = sum( abs(m_std).^chi,2);

end
