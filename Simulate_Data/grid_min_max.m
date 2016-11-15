theta_g_grid = -150:5:100;
theta_c_grid = -250:10:510;
[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);

ds = 1;
grid_lf = rejection_grids_cell{ds,1};


%Vectorize the theta_c, theta_g, and rejection grids,
theta_c_vec = xgrid(:);
theta_g_vec = ygrid(:);

lf_rejection_vec = grid_lf(:);
accepted_values_vec = 1 - lf_rejection_vec;

%Create matrix where each row is an accepted pair of thetas
accepted_thetas = [theta_c_vec( accepted_values_vec == 1) , theta_g_vec( accepted_values_vec == 1) ];

%Find the min and max value of  theta* l
l = [0;1];

max( accepted_thetas * l )
min( accepted_thetas * l )