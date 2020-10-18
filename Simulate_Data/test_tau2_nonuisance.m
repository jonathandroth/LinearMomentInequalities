function [accept, eta, delta, lambda, pval] = test_tau2_nonuisance(tau2, betahat, Sigma, M)
    
    A_nonuisance = [2,0,1];
    A_nonuisance = [A_nonuisance;-A_nonuisance];

    y_nonuisance = A_nonuisance * (betahat - [0;0;1]*tau2) - 3*M;
    X_nonuisance = zeros(size(y_nonuisance));


    SigmaMoments_nonuisance = A_nonuisance * Sigma * A_nonuisance';
    sdMatMoments_nonuisance = diag( sqrt( diag( SigmaMoments_nonuisance ) ) );
    

    y_normalized = sdMatMoments_nonuisance^(-1) * y_nonuisance;
    X_normalized = sdMatMoments_nonuisance^(-1) * X_nonuisance;
    Sigma_normalized = sdMatMoments_nonuisance^(-1) * SigmaMoments_nonuisance * sdMatMoments_nonuisance^(-1);
    
    [reject, eta, delta, lambda, pval] = lp_conditional_test_fn(y_normalized, X_normalized, Sigma_normalized, 0.05);
    
    accept = 1- reject;
end

