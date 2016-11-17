% theta_g_grid = -150:5:100;
% theta_c_grid = -250:10:510;
[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);
%Vectorize the theta_c, theta_g, and rejection grids,
theta_c_vec = xgrid(:);
theta_g_vec = ygrid(:);


confidence_sets_using_grid = NaN(numdatasets,2);
parfor ds = 1:numdatasets
grid_lf = interacted_rejection_grids_cell{ds,1};


grid_lf = grid_lf';
lf_rejection_vec = grid_lf(:);
accepted_values_vec = 1 - lf_rejection_vec;

%Create matrix where each row is an accepted pair of thetas
accepted_thetas = [theta_c_vec( accepted_values_vec == 1) , theta_g_vec( accepted_values_vec == 1) ];

%Find the min and max value of  theta* l
%l = [0;1];
confidence_sets_using_grid(ds,:) = [min( accepted_thetas * l ) , max( accepted_thetas * l )];

end

