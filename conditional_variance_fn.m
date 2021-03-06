

function conditional_variance = conditional_variance_fn(Y, X)

%% Estimate the conditional variance E[ Var(Y_i|Z_i) ] using the method of Abadie Imbens and Zhang (2014)
%Inputs
% Y: an N x K matrix where each row i corresponds with the K-dimensional
% vector Y_i for unit i
% Z: an N x M matrix where each row i correspoinds with the M-dimensional
% vector Z_i for unit
%Outputs:
% conditonal_variance: An estimate of E[ Var(Y_i | Z_i) ]


    T = size(Y,1);
    

    Sigma_X = cov(X);
   
    %If Sigma_X does not have full rank, take an independent subset of the
    %columns of X
    if( rank(Sigma_X) < size(Sigma_X,1) )
        
        [~,licols] = rref( Sigma_X );
        
        X = X(:, licols);
        Sigma_X = Sigma_X(licols, licols);
    end
    
    y_matched = Y( closest_neighbor_indices(  X , Sigma_X), :) ;
    
    conditional_variance = 1/(2*T) * (Y- y_matched)' * (Y - y_matched);

end
