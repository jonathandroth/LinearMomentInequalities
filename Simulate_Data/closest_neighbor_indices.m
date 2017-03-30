%Inputs
%Z_mat: a matrix (nxk), the rows of which represent vectors to be compared
%Sigma: a square k x k covariance matrix

%Outputs
% index_vec: the indices indicating the row that is closest to row i in
% terms of normalized distance
function index_vec = closest_neighbor_indices( Z_mat , Sigma)

numrows = size(Z_mat,1);
index_vec = NaN( numrows,1 );

Sigma_inv = Sigma^(-1);

for row = 1:size(Z_mat,1)
   
    
    Z_t = Z_mat(row,:);
    Z_minus_t = Z_mat( (1:numrows) ~= row, :) ; %all rows of z except for t
    
    
    %Compute the row the minimizes the distance. This row is the row of
    %Z_minus_t
    Z_minus_t_min_row = closest_neighbor_index_helper( Z_t, Z_minus_t, Sigma_inv); 
    
    %Convert the row into a row of Z_mat. In particular, since we omitted
    %the "row"th row, we add 1 if the min row is >= row
    if Z_minus_t_min_row >= row
        Z_minus_t_min_row = Z_minus_t_min_row + 1;
    end
    
    index_vec(row) = Z_minus_t_min_row;
end


end






%Inputs:
%Z_t a row vector (1 x k)
%Z_mat: a (mxk) matrix where each row is a vector to be compared to Z_t
% Sigma_inv: the inverse of a kxk covariance matrix

%Output
% i: the index of the row of z_mat with the shortest normalized distance to
% Z_t

function i = closest_neighbor_index_helper( Z_t, Z_mat, Sigma_inv)

distance_fn = @(Z_s) (Z_t - Z_s) * Sigma_inv * (Z_t - Z_s)' ;
Z_cell = num2cell( Z_mat, 2);
distance = cellfun( distance_fn, Z_cell);
[~,i] = min(distance);

end
