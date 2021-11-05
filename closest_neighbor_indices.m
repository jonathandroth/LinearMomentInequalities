function index_vec = closest_neighbor_indices( Z_mat , Sigma)
% A function for finding the indices of nearest neighbors given a matrix Z
% and covariance matrix Sigma, where distance is in the Mahalanobis
% distance
%Inputs
%Z_mat: a matrix (nxk), the rows of which represent vectors to be compared
%Sigma: a square k x k covariance matrix
%Outputs
% index_vec: the indices indicating the row that is closest to row i in
% terms of normalized distance



%Compute Mahalanobis distance btwn the vectors in Z 
distanceMat = squareform( pdist( Z_mat, 'mahalanobis', Sigma) );

%Replace the diagonol with Inf, so that when we take the min we don't get
%the point itself back
distanceMat( logical(eye(size(distanceMat) ) ) ) = Inf;

[~,index_vec] = min(distanceMat, [] ,2);

end