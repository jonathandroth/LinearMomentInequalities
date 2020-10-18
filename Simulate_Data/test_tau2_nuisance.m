function [accept, eta, delta, lambda, pval] = test_tau2_nuisance(tau2, betahat, Sigma, M)
    A = [1, 1, 0; 0, -2, 1] ;
    A = [A;-A];

    %Create X that treats tau1 as the nuisance parameter
    tauCoefs = [0; 1 ; 0];
    X = A *tauCoefs;
    
    y = A * (betahat - [0;0;1]*tau2) - M;

    SigmaMoments = A * Sigma * A';
    sdMatMoments = diag( sqrt( diag( SigmaMoments ) ) );
    
    y_normalized = sdMatMoments^(-1) * y;
    X_normalized = sdMatMoments^(-1) * X;
    Sigma_normalized = sdMatMoments^(-1) * SigmaMoments * sdMatMoments^(-1);
    
    [reject, eta, delta, lambda, pval] = lp_conditional_test_fn(y_normalized, X_normalized, Sigma_normalized, 0.05);
    
    accept = 1- reject;
end
