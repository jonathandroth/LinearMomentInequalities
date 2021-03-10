%% Add KMS to the path
%addpath('./KMSCode/KMSportable_V3')
addpath(genpath(pwd))
%% Load the data
create_unconditional_moments_1thetac_basic;
%create_unconditional_moments_3thetacs_basic;
%create_unconditional_moments_3thetacs_interacted;

%% Compute identified set bounds using the generated unconditional moments
Y_bar = mean(Y_mat)';
X_c_bar = mean(A_c_combined)';
X_g_bar = mean(A_g_combined)';
X_bar = [X_c_bar, X_g_bar];

mean_g = mean(G_array(1,:,1));
l = [ones(num_F_groups_parameters,1) / num_F_groups_parameters; mean_g];

id_set = cs_linear_delta_lp_fn(Y_bar,X_bar,l,0);


%Compare ID sets using fourmotz and cs_linear_delta_lp_fn
    %Create a change of basis for X so that l corresponds with second column
    M = [null( repmat(l',size(l)) )' ; l'];
    X_tilde_bar = X_bar * M^(-1) ;

    id_set
    [Afm,bfm] = fourmotz(-X_tilde_bar, -Y_bar)
    bfm ./ Afm

 
    
 %% Create unconditional moments at theta_ub when l'delta is the target param
 theta_ub = max(id_set);
 
 %Change of basis so that last column of X_tilde_bar corresponds with l
 M = [null( repmat(l',size(l)) )' ; l'];
 X_tilde_bar = X_bar * M^(-1);
 
 Y_tilde_mat = Y_mat - (theta_ub) * (X_tilde_bar(:,end))';
 X = X_tilde_bar(:,1:(end-1));
 
 %% SImulations using unconditional moments
 alpha = 0.05;
 
 
 numsims = 100;
 nummarkets = 500;
 
 hybrid_test_rejection_vec_unconditional = nan(numsims,1);
 cc_test_rejection_vec_fm_unconditional = nan(numsims,1);
 cc_test_rejection_vec_optimized_unconditional = nan(numsims,1);
 cc_test_rejection_vec_cox_and_shi_code_unconditional = nan(numsims,1);
 cc_test_stat_and_cv_vec_unconditional = nan(numsims,2);
 
% as_test_rejection_vec_unconditional = nan(numsims,1);
 as_test_rejection_vec_unconditional1 = nan(numsims,1);
 as_test_rejection_vec_unconditional2 = nan(numsims,1);
 as_test_rejection_vec_unconditional3 = nan(numsims,1);
 kms_test_rejection_vec_unconditional = nan(numsims,1);
 %Simulations using unconditional moments 
 parfor s = 1:numsims
%for s = 1:numsims
%    s
    rng(s);
    rand_index = randi(size(Y_tilde_mat,1), nummarkets, 1);
    Y_mat_sample = Y_tilde_mat(rand_index, :);
    Sigmahat = cov(Y_mat_sample);
    y_T = mean(Y_mat_sample)'* sqrt(nummarkets);
    D = diag(diag(Sigmahat))^(-0.5);
    y_T = D * y_T;
    X_sim = D* X
    Sigmahat = D*Sigmahat*D;
    
    hybrid_test_rejection_vec_unconditional(s) = lp_hybrid_test_fn( y_T, X_sim, Sigmahat, alpha, alpha/10);
    cc_test_rejection_vec_fm_unconditional(s) = CoxAndShi_FM(y_T,X_sim,Sigmahat,alpha);
    [QLR_stat, ~, CS_crit] = rcc_test_fn(y_T/sqrt(nummarkets), eye(size(y_T,1)) , X_sim, zeros(size(y_T,1),1) , Sigmahat, nummarkets, 1, alpha);
    cc_test_rejection_vec_optimized_unconditional(s) = QLR_stat > CS_crit;
    cc_test_stat_and_cv_vec_unconditional(s,:) = [QLR_stat , CS_crit];
    
    [T_CC,cv_CC, dof_n] = func_subCC(X_sim, -y_T, Sigmahat, alpha);
    cc_test_rejection_vec_cox_and_shi_code_unconditional(s) = (dof_n > 0) * (T_CC > cv_CC);
    
    %Compute AS projection for the second component of X_tilde (l'delta)
    % For AS, we use the full matrix X_tilde_bar, and the oriignal moments
    % that have not been shifted by the ub for l'delta
    Y_mat_sample = Y_mat(rand_index, :);
    Sigmahat = cov(Y_mat_sample);
    y_T = mean(Y_mat_sample)'* sqrt(nummarkets);
    D = diag(diag(Sigmahat))^(-0.5);
    y_T = D * y_T;
    Sigmahat = D*Sigmahat*D;
    X_tilde_T = sqrt(nummarkets)* X_tilde_bar;
    AS_interval1 = projected_AS_or_KMS(y_T, D*X_tilde_T, nummarkets, Sigmahat,[0;1], NaN, 'AS'); %AS interval with standardized moments
    %AS_interval2 = projected_AS_or_KMS(mean(Y_mat_sample)'*sqrt(nummarkets), X_tilde_bar, nummarkets, cov(Y_mat_sample), [0;1], NaN) %AS interval with non-standardized moments
    %AS_interval3 = projected_AS_or_KMS(NaN, X_tilde_bar, NaN, NaN, [0;1], Y_mat_sample); %use the original data rather than sims based off the y_T, X_T
    as_test_rejection_vec_unconditional1(s) = AS_interval1(1) > theta_ub | AS_interval1(2) < theta_ub;
    %as_test_rejection_vec_unconditional2(s) = AS_interval2(1) > theta_ub | AS_interval2(2) < theta_ub;
    %as_test_rejection_vec_unconditional3(s) = AS_interval3(1) > theta_ub | AS_interval3(2) < theta_ub;
    
     kms_interval = projected_AS_or_KMS(y_T, D*X_tilde_T, nummarkets, Sigmahat,[0;1], NaN, 'KMS'); %KMS interval with standardized moments
     kms_test_rejection_vec_unconditional(s) = kms_interval(1) > theta_ub | kms_interval(2) < theta_ub; 
 end    
 
 mean(hybrid_test_rejection_vec_unconditional)
 mean(cc_test_rejection_vec_fm_unconditional)
 mean(cc_test_rejection_vec_optimized_unconditional)
 mean(cc_test_rejection_vec_cox_and_shi_code_unconditional)
% mean(as_test_rejection_vec_unconditional)
 mean(as_test_rejection_vec_unconditional1)
% mean(as_test_rejection_vec_unconditional2)
% mean(as_test_rejection_vec_unconditional3)
 mean(kms_test_rejection_vec_unconditional)
 %% Simulations using conditional moments
 
 hybrid_test_rejection_vec = nan(numsims,1);
 cc_test_rejection_vec_fm = nan(numsims,1);
 cc_test_rejection_vec_optimized = nan(numsims,1);
 cc_test_stat_and_cv_vec = nan(numsims,2);
 cc_test_rejection_vec_cox_and_shi_code= nan(numsims,2);
 
 %Simulations using unconditional moments 
 parfor s = 1:numsims
    rng(s);
    rand_index = randi(size(Y_tilde_mat,1), nummarkets, 1);
    Y_mat_sample = Y_mat(rand_index, :);
    
    A_c_sample = A_c_combined(rand_index, :);
    A_g_sample = A_g_combined(rand_index, :);
    Sigmahat = conditional_variance_fn(Y_mat_sample, [A_c_sample,A_g_sample], 0);
    
    X_bar = [mean(A_c_sample)', mean(A_g_sample)'];
    %Change of basis so that last column of X_tilde_bar corresponds with l
    M = [null( repmat(l',size(l)) )' ; l'];
    X_tilde_bar = X_bar * M^(-1);
 
    Y_bar = mean(Y_mat_sample)';
    Ytilde_bar = Y_bar - (theta_ub) * (X_tilde_bar(:,end));
    X = X_tilde_bar(:,1:(end-1));
    y_T = Ytilde_bar * sqrt(nummarkets);
    D = diag(diag(Sigmahat))^(-0.5);
    y_T = D * y_T;
    X = D * X;
    Sigmahat = D*Sigmahat*D;
    
    hybrid_test_rejection_vec(s) = lp_hybrid_test_fn( y_T, X, Sigmahat, alpha, alpha/10);
    cc_test_rejection_vec_fm(s) = CoxAndShi_FM(y_T,X,Sigmahat,alpha);
    [QLR_stat, ~, CS_crit] = rcc_test_fn(y_T/sqrt(nummarkets), eye(size(y_T,1)) , X, zeros(size(y_T,1),1) , Sigmahat, nummarkets, 1, alpha);
    cc_test_stat_and_cv_vec(s,:) = [QLR_stat , CS_crit];
    cc_test_rejection_vec_optimized(s) = QLR_stat > CS_crit;
    
    [T_CC,cv_CC, dof_n] = func_subCC(X, -y_T, Sigmahat, alpha);
    cc_test_rejection_vec_cox_and_shi_code(s) = (dof_n > 0) * (T_CC > cv_CC);

 end    
 
 mean(hybrid_test_rejection_vec)
 mean(cc_test_rejection_vec_fm)
 mean(cc_test_rejection_vec_optimized)
 mean(cc_test_rejection_vec_cox_and_shi_code)
 

 %% Compare conditional and unconditional
 
 
 rand_index = randi(size(Y_tilde_mat,1), nummarkets, 1);
 Y_mat_sample = Y_mat(rand_index, :);
    
 A_c_sample = A_c_combined(rand_index, :);
 A_g_sample = A_g_combined(rand_index, :);
 Sigmahat_unconditional = cov(Y_mat_sample);
 Sigmahat = conditional_variance_fn(Y_mat_sample, [A_c_sample,A_g_sample], 0);

[ cc_test_stat_and_cv_vec_unconditional(:,1),  cc_test_stat_and_cv_vec(:,1)]
cc_test_stat_and_cv_vec_unconditional(:,1) >  cc_test_stat_and_cv_vec(:,1) 



%% Simulations with unconditional moments and all binding at optimum delta
[eta, deltastar] = test_delta_lp_fn( mean(Y_tilde_mat)', X);
slackness = mean(Y_tilde_mat)' - deltastar * X;

%shift Y_tilde_mat by slackness
Y_tilde_mat_shifted = Y_tilde_mat - slackness';
    %Verify that all moments have ~0 slackness
    [eta, deltastar] = test_delta_lp_fn( mean(Y_tilde_mat_shifted)', X);
    slackness = mean(Y_tilde_mat_shifted)' - deltastar * X;
    max(abs(slackness))
    %Check if FM moments are binding
    [A_fm, b_fm] = findFMMatAndVec(X);
    max(abs(A_fm * mean(Y_tilde_mat_shifted)' - b_fm))

    
 hybrid_test_rejection_vec_unconditional = nan(numsims,1);
 cc_test_rejection_vec_fm_unconditional = nan(numsims,1);
 cc_test_rejection_vec_optimized_unconditional = nan(numsims,1);
 cc_test_stat_and_cv_vec_unconditional = nan(numsims,2);
 cc_test_stat_and_cv_vec_fm_unconditional = nan(numsims,2);
 cc_test_rejection_vec_cox_and_shi_code_unconditional = nan(numsims,1);
 %Simulations using unconditional moments 
for s = 1:100
% parfor s = 1:numsims
    rng(s);
    rand_index = randi(size(Y_tilde_mat_shifted,1), nummarkets, 1);
    Y_mat_sample = Y_tilde_mat_shifted(rand_index, :);
    Sigmahat = cov(Y_mat_sample);
    y_T = mean(Y_mat_sample)'* sqrt(nummarkets);
    D = diag(diag(Sigmahat))^(-0.5);
    y_T = D * y_T;
    X_sim = D*X;
    Sigmahat = D*Sigmahat*D;
    
    hybrid_test_rejection_vec_unconditional(s) = lp_hybrid_test_fn( y_T, X_sim, Sigmahat, alpha, alpha/10);
    [fm_reject, QLR_fm, crit_FM] = CoxAndShi_FM(y_T,X_sim,Sigmahat,alpha);
    cc_test_rejection_vec_fm_unconditional(s) = fm_reject;
    cc_test_stat_and_cv_vec_fm_unconditional(s,:) = [QLR_fm , crit_FM];
    [QLR_stat, ~, CS_crit] = rcc_test_fn(y_T/sqrt(nummarkets), eye(size(y_T,1)) , X_sim, zeros(size(y_T,1),1) , Sigmahat, nummarkets, 1, alpha);
    cc_test_rejection_vec_optimized_unconditional(s) = QLR_stat > CS_crit;
    cc_test_stat_and_cv_vec_unconditional(s,:) = [QLR_stat , CS_crit];
    
    [T_CC,cv_CC, dof_n] = func_subCC(X_sim, -y_T, Sigmahat, alpha);
    cc_test_rejection_vec_cox_and_shi_code_unconditional(s) = (dof_n > 0) * (T_CC > cv_CC);

    
    if(cc_test_rejection_vec_fm_unconditional(s) ~= cc_test_rejection_vec_optimized_unconditional(s))
        warning("FM and non-FM versions don't agree");
        
        
            %Compute AS projection for the second component of X_tilde (l'delta)
    % For AS, we use the full matrix X_tilde_bar, and the oriignal moments
    % that have not been shifted by the ub for l'delta
    Y_mat_sample = Y_mat(rand_index, :) - slackness';
    Sigmahat = cov(Y_mat_sample);
    y_T = mean(Y_mat_sample)'* sqrt(nummarkets);
    D = diag(diag(Sigmahat))^(-0.5);
    y_T = D * y_T;
    Sigmahat = D*Sigmahat*D;
    AS_interval1 = projected_AS_or_KMS(y_T, D*X_tilde_bar, nummarkets, Sigmahat,[0;1], NaN); %AS interval with standardized moments
    as_test_rejection_vec_unconditional1(s) = AS_interval1(1) > theta_ub | AS_interval1(2) < theta_ub;

    end
 end    
 
 mean(hybrid_test_rejection_vec_unconditional)
 mean(cc_test_rejection_vec_fm_unconditional)
 mean(cc_test_rejection_vec_optimized_unconditional)
 mean(cc_test_rejection_vec_cox_and_shi_code_unconditional)
 
 mean(as_test_rejection_vec_unconditional1)

 %% Simulations using shifted conditional moments
 
 %% Simulations using conditional moments
 
 hybrid_test_rejection_vec = nan(numsims,1);
 cc_test_rejection_vec_fm = nan(numsims,1);
 cc_test_rejection_vec_optimized = nan(numsims,1);
 cc_test_stat_and_cv_vec = nan(numsims,2);
 
 %Simulations using unconditional moments 
 parfor s = 1:numsims
    rng(s);
    rand_index = randi(size(Y_tilde_mat,1), nummarkets, 1);
    Y_tilde_mat_sample = Y_tilde_mat_shifted(rand_index, :);
    
    A_c_sample = A_c_combined(rand_index, :);
    A_g_sample = A_g_combined(rand_index, :);
    Sigmahat = conditional_variance_fn(Y_tilde_mat_sample, [A_c_sample,A_g_sample], 0);
    
    X_bar = [mean(A_c_sample)', mean(A_g_sample)'];
    %Change of basis so that last column of X_tilde_bar corresponds with l
    M = [null( repmat(l',size(l)) )' ; l'];
    X_tilde_bar = X_bar * M^(-1);
 
    %Y_bar = mean(Y_mat_sample)';
    %Ytilde_bar = Y_bar - (theta_ub) * (X_tilde_bar(:,end));
    Ytilde_bar = mean(Y_tilde_mat_sample)';
    X = X_tilde_bar(:,1:(end-1));
    y_T = Ytilde_bar * sqrt(nummarkets);
    D = diag(diag(Sigmahat))^(-0.5);
    y_T = D * y_T;
    X = D * X;
    Sigmahat = D*Sigmahat*D;
    
    hybrid_test_rejection_vec(s) = lp_hybrid_test_fn( y_T, X, Sigmahat, alpha, alpha/10);
    cc_test_rejection_vec_fm(s) = CoxAndShi_FM(y_T,X,Sigmahat,alpha);
    [QLR_stat, ~, CS_crit] = rcc_test_fn(y_T/sqrt(nummarkets), eye(size(y_T,1)) , X, zeros(size(y_T,1),1) , Sigmahat, nummarkets, 1, alpha);
    cc_test_stat_and_cv_vec(s,:) = [QLR_stat , CS_crit];
    cc_test_rejection_vec_optimized(s) = QLR_stat > CS_crit;
 end    
 
 mean(hybrid_test_rejection_vec)
 mean(cc_test_rejection_vec_fm)
 mean(cc_test_rejection_vec_optimized)