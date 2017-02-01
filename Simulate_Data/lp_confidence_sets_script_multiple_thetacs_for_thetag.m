%This partial script calculates confidence sets for theta_g

%Assumes that the confidence sets for theta_c have already been calculated,
%so does not need to load the moment function

l = [zeros( size(F_group_cell,1) ,1) ; 1];
xlim_graph = [-150;100];
filename_graph =  'Theta_g Rejection Probabilities';
xlabel_graph = 'Theta_g';
data_output_dir = strcat(data_output_dir, 'Theta_g/');

lp_compute_confidence_sets;
lp_graphs;