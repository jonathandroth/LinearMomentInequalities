%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Operations on the data for inference
%   For the sample problems, we only need to use mbar_std
%   For the MR problem, we only need mData and xi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mbar_std,mData,xi] = dataOperations_function(theta,data,kappa)
% theta denotes the parameter of interest.
% data denotes the dataset.
% kappa denotes the tuning parameter.

% sample size;
n = size(data,1);

% dimension of X
dimSX = size(data,2)/2;

% computes auxiliary information from data;
dataP11 = data(:,1:dimSX); % represents P(A_1=1,A_2=1)
dataP10 = data(:,dimSX+1:2*dimSX); % represents P(A_1=1,A_2=0)

% creates auxiliary parameters
auxParameter = [0,theta(3:end)]; % defines beta_q for q=1,...,d_X
auxParameter1 = theta(1) - auxParameter; % defines theta_1 - beta_q for q=1,...,d_X
auxParameter2 = theta(2) - auxParameter; % defines theta_2 - beta_q for q=1,...,d_X

% computes observations whose expectation should be Eq. (5.1)
mData_1q  = dataP11 - repmat( (1-auxParameter1).*(1-auxParameter2) ,n,1); % Equalities in m_{1,q}
mData_2q = dataP10 - repmat( auxParameter2.*(1-auxParameter1),n,1); % Inequalities in m_{2,q}
mData_3q = repmat( auxParameter2 ,n,1) - dataP10; % Inequalities in m_{3,q}
mData = [mData_2q, mData_3q, mData_1q]; % defines moment (in)equalities data. Note: inequalities should appear first

% compute studentized sample averages of mData 
epsilon = 0.000001; % introduces this parameter to avoid division by zero
mbar_std  = mean(mData)./(std(mData) + epsilon);

% Additional parameter needed in DR, PR, and MR inference
xi = (1/kappa)*sqrt(n)*mbar_std; % Slackness measure in GMS

end