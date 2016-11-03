

function conditional_variance = conditional_variance_fn(y, X, diagonol)


    T = size(y,1);
    

    Sigma_X = cov(X);
    
    if(diagonol)
        Sigma_X = diag( diag( Sigma_X) );
        %display('Using Diagonal Matrix');
    
    elseif( rank(Sigma_X) < size(Sigma_X,1) )
        
        Sigma_X = Sigma_X + 0.01 * diag( diag( Sigma_X) );
    end
    
    y_matched = y( closest_neighbor_indices(  X , Sigma_X), :) ;
    
    conditional_variance = 1/(2*T) * (y- y_matched)' * (y - y_matched);

end