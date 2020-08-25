

%betahat = [-1;1;3];
%betahat = [-1;1;-2]; %use this for cases where nuisance and no-nuisance not analytically the same

betahat = [0;1;4];
%betahat = [-8;1;6];
%betahat = [-20;1;50];

%Sigma = eye(length(betahat));
%Sigma = [1,0.5,0.5; 0.5, 1, 0.5; 0.5, 0.5, 1];
Sigma = [2,0.5,0.5; 0.5, 1, 0.5; 0.5, 0.5, 1];

M = 0.5;
tau2 = 0;


taugrid = -4:.2:4;

tauCI_nonuisance = NaN(size(taugrid));
tauCI_nuisance = NaN(size(taugrid));

eta_nonuisance = NaN(size(taugrid));
eta_nuisance = NaN(size(taugrid));

pval_nonuisance = NaN(size(taugrid));
pval_nuisance = NaN(size(taugrid));


for(i = 1:length(taugrid) )
    [tauCI_nonuisance(i), eta_nonuisance(i), ~ , ~ , pval_nonuisance(i)] = test_tau2_nonuisance(taugrid(i), betahat, Sigma, M);
    [tauCI_nuisance(i), eta_nuisance(i), ~, ~, pval_nuisance(i)] = test_tau2_nuisance(taugrid(i), betahat, Sigma, M);

end
%test_tau2_nuisance(tau2, betahat, Sigma, M)
%test_tau2_nonuisance(tau2, betahat, Sigma, M)

tauCI_nonuisance
tauCI_nuisance


%% To do
% You should try comparing the ratios of etas under the two programs, with
% non-homogenous sigmas

%Does the ratio of etas line up with the factor that the theory predicts in
%the lower bound?
