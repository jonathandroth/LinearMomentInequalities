

function conditional_variance = conditional_variance_fn(Y, X, diagonol)

%% Estimate the conditional variance E[ Var(Y_i|Z_i) ] using the method of Abadie Imbens and Zhang (2014)
%Inputs
% Y: an N x K matrix where each row i corresponds with the k-dimensional
% vector Y_i for unit i
% Z: an N x M matrix where each row i correspoinds with the m-dimensional
% vector Z_i for unit
% Sigma: the (scaled) estimate of E[ Var(Y_i|Z_i) ]
% l: The target parameter is l'delta
% alpha: the size of the test (e.g. 0.05 for 5% significance)
% lf_cv: The least favorable CV outputted by lf_critical_value_fn
%Outputs:
% cs: a two-dimensional vector with the lower and upper bounds of the
% confidence set
% slack_lb: the slackness of the moments at the lb bound of the cs
% slack_ub: the slackness of the moments at the ub bound of the cs
%Notes
% If y and X are sample averages, Sigma should be E[ Var(Y|X) ] /N
% If Sigma is E[ Var(Y|X) ], then y and X should be averages scaled by sqrt(N)
% Sigma used here should be the same as the sigma used for lf_cv



    T = size(Y,1);
    

    Sigma_X = cov(X);
    
    if(diagonol)
        Sigma_X = diag( diag( Sigma_X) );
        %display('Using Diagonal Matrix');
    
    %If Sigma_X does not have full rank, take an independent subset of the
    %columns of X
    elseif( rank(Sigma_X) < size(Sigma_X,1) )
        
        [~,licols] = rref( Sigma_X );
        
        X = X(:, licols);
        Sigma_X = Sigma_X(licols, licols);
    end
    
    y_matched = Y( closest_neighbor_indices(  X , Sigma_X), :) ;
    
    conditional_variance = 1/(2*T) * (Y- y_matched)' * (Y - y_matched);

end
