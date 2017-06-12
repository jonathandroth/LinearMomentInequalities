
%% Specify parameters
lp_set_parameters

if( graphs_only == 0 )

%% Import data and calculate moments and covariance matrices
lp_create_moments_and_covariances
%% Calculate the confidence sets and the identified set
lp_compute_confidence_sets

end
%% Create graphs
lp_graphs
