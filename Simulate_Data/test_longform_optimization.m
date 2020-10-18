%% 2 period example
numPre = 1;
numPost = 1;

A_SD = [ 1, 1; -1,-1];
M = 1;
K = 1;
taubarK = 0;
alpha = 0.05;

betahat = [1; -5]; 
%Sigma_betahat = [1,.5;1,.5] ;
Sigma_betahat = [2,.5;1,.5] ;

[X,Y, Sigma] = create_X_Y_Sigma_for_longform(numPre, numPost, A_SD, M, K, taubarK, betahat, Sigma_betahat);

[eta_star, delta_star, lambda, error_flag]  = test_delta_lp_fn( Y, X, Sigma)
%lp_conditional_test_fn( Y, X, Sigma, alpha)

%% 3 period example
numPre = 1;
numPost = 2;

A_SD = [ 1, 1,0 ; -1,-1, 0;
         0,-2,1; 0, 2, -1];
M = 1;
K = 2;
taubarK = 0;
alpha = 0.05;

betahat = [0; 1; 20]; 

Sigma_betahat = diag([2;1;1]) ;

[X,Y, Sigma] = create_X_Y_Sigma_for_longform(numPre, numPost, A_SD, M, K, taubarK, betahat, Sigma_betahat);

[eta_star, delta_star, lambda, error_flag]  = test_delta_lp_fn( Y, X, Sigma)

SigmaDiag = diag( diag(Sigma) )
SigmaDiag = SigmaDiag +  diag( diag(SigmaDiag) == 0 );
Y_normalized = SigmaDiag^(-1/2) * Y;
X_normalized = SigmaDiag^(-1/2) * X;
Sigma_normalized = SigmaDiag^(-1/2)  * Sigma * SigmaDiag^(-1/2) ;

lp_conditional_test_fn( Y_normalized, X_normalized, Sigma_normalized, alpha)
%you're having a problem in the dual approach; possibly has to do with the
%fact that you have zero values for the sigmas. You also changed the
%test_delta_lp and lp_conditional_test functions to no longer assume sigma
%is a correlation matrix, and that may have messed things up too
%% function def
function [X,Y, Sigma] = create_X_Y_Sigma_for_longform(numPre, numPost, A, M, K, taubarK, betahat, Sigma_betahat)

Y = [betahat; -betahat; -M * ones(size(A,1),1) ; taubarK; taubarK];

X_betahat_moments_positive_direction = [ eye(numPre), zeros(numPre, 2*numPost); ...
                                         zeros(numPost, numPre), eye(numPost), eye(numPost)];
                                     
X_betahat_moments_negative_direction = -X_betahat_moments_positive_direction;

X_Adelta = [ A , zeros( size(A,1), numPost) ];

e_k = double( 1:numPost == K );
X_tauK = [ zeros(1, numPre), zeros(1, numPost), e_k ];

X = [X_betahat_moments_positive_direction; X_betahat_moments_negative_direction; X_Adelta; X_tauK; -X_tauK];

l = [ eye(numPre + numPost); -eye(numPre + numPost); zeros( size(X,1) - 2* (numPre+ numPost) , numPre + numPost)];
Sigma = l * Sigma_betahat * l';
%Sigma_Betahats = [ eye(numPre + numPost); -eye(numPre + numPost)] * Sigma_betahat * [ eye(numPre + numPost); -eye(numPre + numPost)]' ;

%Sigma = [ Sigma_Betahats, zeros( numPre + numPost,  length(Y) - 2*(numPre + numPost)); ...
%          zeros( length(Y) - 2*(numPre + numPost)   , length(Y)  ) ];
end