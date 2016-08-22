function conditional_variance = conditional_variance_fn(y, X)


    T = size(y,1);
    
    Sigma_X = cov(X);
    y_matched = y( closest_neighbor_indices(  X , Sigma_X), :) ;
    
    conditional_variance = 1/(2*T) * (y- y_matched)' * (y - y_matched);

end