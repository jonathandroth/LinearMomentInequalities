%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Defines the sample criterion function according to MMM in Eq. (4.2);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = Qn_function(theta_to_min,theta_H0,coordinate,data,p)
% theta_to_min denotes the subvector of theta that should be minimized;
% theta_H0 denotes the subvector of theta that is fixed according to H0;
% coodinate denotes the coordinate of interest (1 or 2);
% data denotes the dataset.
% p denotes the number of moment inequalities (which should appear first);

% defines the parameter vector;
if coordinate ==1
    theta = [theta_H0,theta_to_min];
else % coordinate = 2;
    theta = [theta_to_min(1),theta_H0,theta_to_min(2:end)];
end

% determines sample size;
n = size(data,1);

% studentizes and averages the data;
[mbar_std,~,~] = dataOperations_function(theta,data,1); % No use for kappa, set to one.

% computes the sample criterion function;
value = S_function(sqrt(n)*mbar_std,p);