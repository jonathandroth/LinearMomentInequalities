%% Load the data
%XX Load the saves unconditional moments
% You also need to either save l or load the G_array that it's created from
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
 
 Y_tilde_mat = Y_mat - (theta_ub+1) * (X_tilde_bar(:,end))';
 X = X_tilde_bar(:,1:(end-1));
 
 %% SImulations
 alpha = 0.05;
 
 
 numsims = 100;
 nummarkets = 500;
 
 hybrid_test_rejection_vec = nan(numsims,1);
 cc_test_rejection_vec_fm = nan(numsims,1);
  cc_test_rejection_vec_optimized = nan(numsims,1);
 
 parfor s = 1:numsims
    rng(s);
    rand_index = randi(size(Y_tilde_mat,1), nummarkets, 1);
    Y_mat_sample = Y_tilde_mat(rand_index, :);
    Sigmahat = cov(Y_mat_sample);
    y_T = mean(Y_mat_sample)'* sqrt(nummarkets);
    D = diag(diag(Sigmahat))^(-0.5);
    y_T = D * y_T;
    Sigmahat = D*Sigmahat*D;
    
    hybrid_test_rejection_vec(s) = lp_hybrid_test_fn( y_T, X, Sigmahat, alpha, alpha/10);
    cc_test_rejection_vec_fm(s) = CoxAndShi_FM(y_T,X,Sigmahat,alpha);
    [QLR_stat, ~, CS_crit] = rcc_test_fn(y_T/sqrt(nummarkets), eye(size(y_T,1)) , X, zeros(size(y_T,1),1) , Sigmahat, nummarkets, 1, alpha);
    cc_test_rejection_vec_optimized(s) = QLR_stat > CS_crit;
 end    
 
 mean(hybrid_test_rejection_vec)
 mean(cc_test_rejection_vec_fm)
 mean(cc_test_rejection_vec_optimized)
