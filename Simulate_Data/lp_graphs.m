%This script creates graphs for the LP confidence sets runs

%% Create graphs
load(strcat( data_output_dir, dirname, 'Interacted_Moments/confidence_sets_lp'));
load(strcat( data_output_dir, dirname, 'Interacted_Moments/identified_set_bounds'));

gridpoints = 1000;
l_theta_grid = linspace(min(confidence_sets_using_c_alpha(:,1)) -1 ,...
                        max(confidence_sets_using_c_alpha(:,2)) + 1, gridpoints );
                    
rejection_grid_c_alpha = NaN(gridpoints,1);
rejection_grid_c_lp_alpha = NaN(gridpoints,1);

for i = 1:gridpoints
    
    l_theta = l_theta_grid(i);
    
    rejection_grid_c_alpha(i,1) = 1 -mean( (confidence_sets_using_c_alpha(:,1) <= l_theta) & ... 
                                   (l_theta <= confidence_sets_using_c_alpha(:,2)) );
    rejection_grid_c_lp_alpha(i,1) = 1 -mean( (confidence_sets_using_c_lp_alpha(:,1) <= l_theta) & ... 
                                   (l_theta <= confidence_sets_using_c_lp_alpha(:,2)) );
    
                               
end


plot( repmat( l_theta_grid', 1, 2), [rejection_grid_c_alpha, rejection_grid_c_lp_alpha ]);
line( [identified_set_bounds(1);identified_set_bounds(1)] ,[0;1], 'LineStyle', '--', 'Color',  'r');
line( [identified_set_bounds(2);identified_set_bounds(2)] ,[0;1], 'LineStyle', '--', 'Color',  'r');


legend( 'LF','LF (modified)', 'Identified Set Boundary',  'Location','eastoutside' );
ylabel('Rejection Probability');


%If manual bounds are specified for the x-axis limit, impose these
if( exist('xlim_graph') ==1)
    xlim( xlim_graph )    
end

%If no xlabel specified, do 'l * theta*
if( exist('xlabel_graph') == 0)
    xlabel_graph = 'l * theta';    
end

xlabel(xlabel_graph);

%If filename not specified, assume it's means
if( exist('filename_graph') ==0)
       filename_graph =  'Mean Weight Rejection Probabilities';
end


saveas( gcf, strcat(figures_output_dir,filename_graph ), 'epsc');


%options = optimoptions('linprog', 'Display', 'final' );
%linprog( 1, [1;-1], [2,-3],[],[],[],[], [],options )

display('Script complete');