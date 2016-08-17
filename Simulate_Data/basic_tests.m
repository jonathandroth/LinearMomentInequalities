% This function does the LF, RSW, conditional, and hybrid tests given:

% g_T: the vector of moments
% Sigma: the covariance mat
% Z_draws: a matrix of normal draws to do the simulations
% alpha, beta: the size of the tests

function [test_lf, test_rsw, test_conditional, test_hybrid] = basic_tests(g_T, Sigma, Z_draws, alpha, beta)  

        %Calculate the normalized maximum (used for multiple test)
        R_gt = R_gt_sigma( g_T, Sigma);
        
        %Do the LF test
        cutoff_lf = c_lf(Sigma, alpha, Z_draws);
        test_lf =  ( R_gt > cutoff_lf );
        
        %Do the RSW test
        [cutoff_RSW , mu_tilde] = c_RSW( g_T, Sigma, alpha, beta, Z_draws);
        test_rsw = ( R_gt > cutoff_RSW ) .* max( mu_tilde == 0) ;
        
        %Do the conditional test
        
        %cutoff_conditional = c_conditional(g_T, Sigma, alpha);
        %test_conditional = (R_gt > cutoff_conditional);
        
        [~, T_conditional_integrated] = c_conditional(g_T, Sigma, alpha);
        %test_conditional = T_conditional < alpha;
        test_conditional = T_conditional_integrated < alpha;
        
        %Do the hybrid test
        cutoff_lf_beta = c_lf(Sigma, beta, Z_draws);
        test_hybrid = max( T_conditional_integrated < (alpha - beta) , R_gt > cutoff_lf_beta );
        
end
        