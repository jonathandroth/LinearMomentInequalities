%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Defines the sample criterion function for DR or PR inference
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = Qn_MR_function(theta_to_min,theta_H0,coordinate,data,kappa,p,k,W2_AA,MR_type)
% theta_to_min denotes the subvector of theta that should be minimized.
% theta_H0 denotes the subvector of theta that is fixed according to H0.
% coodinate denotes the coordinate of interest (1 or 2).
% data denotes the dataset.
% kappa denotes the tuning parameter using by DR or PR.
% p denotes the number of moment inequalities (which should appear first).
% k denotes the number of moment (in)equalities.
% W2_AA is a vector of random variables used to implement the (multiplier) bootstrap.
% MR_type indicates the type of resampling, i.e., DR or PR.

if coordinate == 1
    theta = [theta_H0,theta_to_min];
else % coodinate = 2;
    theta = [theta_to_min(1),theta_H0,theta_to_min(2:end)];
end

n = size(data,1);
[~,mData,xi] = dataOperations_function(theta,data,kappa);

if MR_type == 1 % DR method;
    value = S_function(W2_AA*zscore(mData,1)/sqrt(n) + repmat(phi_function(xi,p,k-p),size(W2_AA,1),1),p);
else % PR method;
    value = S_function(W2_AA*zscore(mData,1)/sqrt(n) + repmat(xi,size(W2_AA,1),1),p);
end;