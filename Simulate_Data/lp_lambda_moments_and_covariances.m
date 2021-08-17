%% Import data and calculate moments and covariance matrices

%lambda_vec = (0.01:10:20.01)';
%lambda_vec = (0.00:0.01:0.01)';
%lambda_vec = (0.01:10:10.01)';
%lambda_vec = (0.01:.2:5.01)';



conditional_test = NaN(numdatasets, size(lambda_vec,1));
hybrid_test = NaN(size(conditional_test));
lf_test_original = NaN(size(conditional_test));
lf_test_modified = NaN(size(conditional_test));

cc_test = NaN(size(conditional_test));
rcc_test = NaN(size(conditional_test));


timing_vec_lf_lambda = zeros(size(conditional_test));
timing_vec_lfp_lambda = zeros(size(conditional_test));
timing_vec_conditional_lambda = zeros(size(conditional_test));
timing_vec_hybrid_lambda = zeros(size(conditional_test));
timing_vec_rcc_lambda = zeros(size(conditional_test));
timing_vec_cc_lambda = zeros(size(conditional_test));


lambda_identified_set = NaN( size(lambda_vec) );
lambda_identified_set_zerocutoff = NaN( size(lambda_vec) );

i = 1;


for lambda = lambda_vec'
%lambda_true = .386

lambda

lp_create_moments_and_covariances


%The create_moments_and_covariances script creates a new timing vec for
%each iteration (e.g. timing_vec_lf) that records the time to compute the
%conditioanl covariance for each ds. We save these

timing_vec_lf_lambda(:,i) = timing_vec_lf; 
timing_vec_lfp_lambda(:,i) = timing_vec_lfp;
timing_vec_conditional_lambda(:,i) = timing_vec_conditional;
timing_vec_hybrid_lambda(:,i) = timing_vec_hybrid;
timing_vec_cc_lambda(:,i) = timing_vec_cc;
timing_vec_rcc_lambda(:,i) = timing_vec_rcc;



nummoments = size( y_T_cell{ds,1} , 1);
rng(0);
Z_draws_interacted = randn(nummoments, 10000);
%CHANGE THIS BACK WHEN DONE DEBUGGING
%for ds = 1:numdatasets
parfor ds = 1:numdatasets
X_T = X_T_cell{ds};
y_T = y_T_cell{ds};
Sigma = Sigma_conditional_cell{ds};

%Do LF-based tests
ticLFPCV = tic;
c_alpha = c_lf(Sigma, alpha, Z_draws_interacted);
runTimeLFPCV = toc(ticLFPCV);

ticLFCV = tic;
[c_lp_alpha, eta_draws] = c_lf_lp(X_T,Z_draws_interacted(:,1:numsims_lp),Sigma,alpha);
runTimeLFCV = toc(ticLFCV);

ticTestStat = tic;
eta = test_delta_lp_fn( y_T, X_T);
runTimeTestStat = toc(ticTestStat); 

lf_test_original(ds,i) = eta > c_alpha;
lf_test_modified(ds,i) = eta > c_lp_alpha;

timing_vec_lfp_lambda(ds,i) = timing_vec_lfp_lambda(ds,i) + runTimeLFPCV + runTimeTestStat;
timing_vec_lf_lambda(ds,i) = timing_vec_lf_lambda(ds,i) + runTimeLFCV + runTimeTestStat;

ticConditional = tic;
conditional_test(ds,i) = lp_conditional_test_fn( y_T, X_T, Sigma, alpha);
runTimeConditional = toc(ticConditional);

ticHybrid = tic;
hybrid_test(ds,i) = lp_hybrid_test_fn( y_T, X_T, Sigma, alpha, alpha/10, eta_draws);
runTimeHybrid = toc(ticHybrid) + runTimeLFCV; %hybrid uses the draws from the lf cv


timing_vec_conditional_lambda(ds,i) = timing_vec_conditional_lambda(ds,i) + runTimeConditional;
timing_vec_hybrid_lambda(ds,i) = timing_vec_hybrid_lambda(ds,i) + runTimeHybrid;

%Run RCC test for specs with <10 params or <100 moments
%Otherwise, use CC test and give upper bound for RCC based on best
%refinement
if size(X_T,2) < 10 || size(X_T,1) < 100
  try
  [T_RCC,cv_RCC,cv_CC,~,~, dof_n, toc_RCC, toc_CC] = func_subRCC(X_T, -y_T, Sigma, alpha);
  cc_test(ds,i) = (dof_n > 0) * (T_RCC > cv_CC);
  rcc_test(ds,i) = (dof_n > 0) * (T_RCC > cv_RCC);
  
  timing_vec_cc_lambda(ds,i) = timing_vec_cc_lambda(ds,i) + toc_CC;
  timing_vec_rcc_lambda(ds,i) = timing_vec_rcc_lambda(ds,i) + toc_RCC;
  
  catch
  [T_CC,cv_CC, dof_n, toc_CC] = func_subCC(X_T, -y_T, Sigma, alpha);
  cc_test(ds,i) = (dof_n > 0) * (T_CC > cv_CC);
  rcc_test(ds,i) = ((dof_n > 0) * (T_CC > cv_CC)) | ( (dof_n ==1) & (T_CC > chi2inv(1-2*alpha,dof_n)) ); %right now set RCC to 1 if dof_n == 1 and T_CC is above the 1-2*alpha cv  
  
  timing_vec_cc_lambda(ds,i) = timing_vec_cc_lambda(ds,i) + toc_CC;
  timing_vec_rcc_lambda(ds,i) = timing_vec_rcc_lambda(ds,i) + toc_CC;

  end
else
    [T_CC,cv_CC, dof_n, toc_CC] = func_subCC(X_T, -y_T, Sigma, alpha);
    cc_test(ds,i) = (dof_n > 0) * (T_CC > cv_CC);
    rcc_test(ds,i) = ((dof_n > 0) * (T_CC > cv_CC)) | ( (dof_n ==1) & (T_CC > chi2inv(1-2*alpha,dof_n)) ); %right now set RCC to 1 if dof_n == 1 and T_CC is above the 1-2*alpha cv
    
    timing_vec_cc_lambda(ds,i) = timing_vec_cc_lambda(ds,i) + toc_CC;
    timing_vec_rcc_lambda(ds,i) = timing_vec_rcc_lambda(ds,i) + toc_CC;

end



end

i = i+1;
end

%% Save the results

    ds_dir = strcat( data_output_dir, 'Interacted_Moments/');
    mkdir(ds_dir);
    save( strcat(ds_dir, 'lambda_results'), ...
        'conditional_test', 'lf_test_original', 'lf_test_modified', ...
        'hybrid_test', 'lambda_vec', 'cc_test', 'rcc_test',...
        'timing_vec_lf_lambda', 'timing_vec_lfp_lambda',...
        'timing_vec_conditional_lambda', 'timing_vec_hybrid_lambda',...
        'timing_vec_cc_lambda', 'timing_vec_rcc_lambda');