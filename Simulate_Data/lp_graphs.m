 %This script creates graphs for the LP confidence sets runs

%% Create graphs
load(strcat( data_output_dir, dirname, 'Interacted_Moments/confidence_sets_lp'));
load(strcat( data_output_dir, dirname, 'Interacted_Moments/identified_set_bounds_zerocutoff'));
load(strcat( data_output_dir, dirname, 'Interacted_Moments/identified_set_bounds'));


%load(strcat( data_output_dir, dirname, 'Interacted_Moments/identified_set_bounds_zerocutoff'));
%identified_set_bounds = identified_set_bounds_zerocutoff;

addpath('./breakxaxis')
addpath('./export-fig')
%gridpoints = 5000;
%l_theta_grid = linspace(min(confidence_sets_using_c_alpha(:,1)) -1 ,...
%                        max(confidence_sets_using_c_alpha(:,2)) + 1, gridpoints );

gridpoints = length(beta0_grid);
l_theta_grid = beta0_grid;

rejection_grid_c_alpha = NaN(gridpoints,1);
rejection_grid_c_lp_alpha = NaN(gridpoints,1);

for i = 1:gridpoints
    
    l_theta = l_theta_grid(i);
    
    rejection_grid_c_alpha(i,1) = 1 -mean( (confidence_sets_using_c_alpha(:,1) <= l_theta) & ... 
                                   (l_theta <= confidence_sets_using_c_alpha(:,2)) );
    rejection_grid_c_lp_alpha(i,1) = 1 -mean( (confidence_sets_using_c_lp_alpha(:,1) <= l_theta) & ... 
                                   (l_theta <= confidence_sets_using_c_lp_alpha(:,2)) );
    
                               
end

getrow = @(x,row) x(row,:);

hold('on');

%p = plot( repmat( l_theta_grid', 1, 2), [rejection_grid_c_alpha, rejection_grid_c_lp_alpha ]);

l1 = plot(l_theta_grid', rejection_grid_c_alpha, 'Color', getrow( get(gca,'colororder'),1 ) )
l2 = plot(l_theta_grid', rejection_grid_c_lp_alpha, 'Color', getrow( get(gca,'colororder'),2 ) )
l3 = plot( beta0_grid , rejection_grid_conditional, 'Color', getrow( get(gca,'colororder'),4 ), 'LineStyle', ':')
l4 = plot( beta0_grid , rejection_grid_hybrid, 'Color', getrow( get(gca,'colororder'),3 ), 'LineStyle', '-.' )

l5 = plot( [identified_set_bounds(1);identified_set_bounds(1)] ,[0;1], 'LineStyle', '--', 'Color',  'r');
l6 = plot( [identified_set_bounds(2);identified_set_bounds(2)] ,[0;1], 'LineStyle', '--', 'Color',  'r');

%p.Color(4) = 0.7;
% l1.Color(4) = 0.7;
% l2.Color(4) = 0.7;
% l3.Color(4) = 0.7;
% l4.Color(4) = 0.7;

legend( 'LFP','LF', 'Conditional', 'Hybrid', 'Identified Set Boundary',  'Location',....
    'southoutside', 'Orientation', 'horizontal' );
ylabel('Rejection Probability');

%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 2.5]);  

% %If no xlabel specified, do 'l * theta*
% if( exist('xlabel_graph') == 0)
%     xlabel_graph = 'l * theta';    
% end
% 
% xlabel(xlabel_graph);


if( exist('xlim_graph') ==1)
    %If manual bounds are specified for the x-axis limit, impose these
    xlim( xlim_graph )    
    
    %Impose tick width if specified; otherwise 5
    if(exist('xtick_width'))
        set(gca,'XTick',[xlim_graph(1):xtick_width:xlim_graph(2)])
    end

end

    
%Create a break in the x-axis if specified
    if(exist('xsplit_graph'))
        breakxaxis(xsplit_graph)
    end

    


%If filename not specified, assume it's means
if( exist('filename_graph') ==0)
       filename_graph =  'Mean_Weight_Rejection_Probabilities';
end

set(findall(gcf,'-property','FontSize'),'FontSize',14);

set(findall(gcf, 'Type', 'Line'),'LineWidth',3); %Linewidth for plot lines

export_fig(strcat(figures_output_dir,filename_graph,'.pdf'));

clf

%options = optimoptions('linprog', 'Display', 'final' );
%linprog( 1, [1;-1], [2,-3],[],[],[],[], [],options )

%% Create plot comparing sCC and rCC test of Cox & Shi to hybrid

hold('on');

%p = plot( repmat( l_theta_grid', 1, 2), [rejection_grid_c_alpha, rejection_grid_c_lp_alpha ]);

l1 = plot( beta0_grid , rejection_grid_hybrid, 'Color', getrow( get(gca,'colororder'),3 ), 'LineStyle', '-.' )
l2 = plot( beta0_grid , rejection_grid_cc, 'Color', getrow( get(gca,'colororder'),6 ), 'LineStyle', '-' )
l3 = plot( beta0_grid , rejection_grid_rcc, 'Color', getrow( get(gca,'colororder'),7 ), 'LineStyle', ':' )

l4 = plot( [identified_set_bounds(1);identified_set_bounds(1)] ,[0;1], 'LineStyle', '--', 'Color',  'r');
l5 = plot( [identified_set_bounds(2);identified_set_bounds(2)] ,[0;1], 'LineStyle', '--', 'Color',  'r');

legend( 'Hybrid', 'sCC', 'sRCC', 'Identified Set Boundary',  'Location',....
    'southoutside', 'Orientation', 'horizontal' );
ylabel('Rejection Probability');

%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 2.5]);  

% %If no xlabel specified, do 'l * theta*
% if( exist('xlabel_graph') == 0)
%     xlabel_graph = 'l * theta';    
% end
% 
% xlabel(xlabel_graph);


if( exist('xlim_graph') ==1)
    %If manual bounds are specified for the x-axis limit, impose these
    xlim( xlim_graph )    
    
    %Impose tick width if specified; otherwise 5
    if(exist('xtick_width'))
        set(gca,'XTick',[xlim_graph(1):xtick_width:xlim_graph(2)])
    end

end

    
%Create a break in the x-axis if specified
    if(exist('xsplit_graph'))
        breakxaxis(xsplit_graph)
    end

    


%If filename not specified, assume it's means
if( exist('filename_graph') ==0)
       filename_graph =  'Mean_Weight_Rejection_Probabilities';
end

set(findall(gcf,'-property','FontSize'),'FontSize',14);

set(findall(gcf, 'Type', 'Line'),'LineWidth',3); %Linewidth for plot lines

export_fig(strcat(figures_output_dir,filename_graph, '_compare_to_Cox_and_Shi','.pdf'));

clf


%% Create rows for table that shows where we achieve various power thresholds

[lb_row_05,ub_row_05] = ub_and_lb_row_for_table( beta0_grid, l_theta_grid,...
                                                     rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,...
                                                     identified_set_bounds,...
                                                     0.05);
                                                 
 [lb_row_50,ub_row_50] = ub_and_lb_row_for_table( beta0_grid, l_theta_grid,...
                                                     rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,...
                                                     identified_set_bounds,...
                                                     0.50);
                                                 
[lb_row_95,ub_row_95] = ub_and_lb_row_for_table( beta0_grid, l_theta_grid,...
                                                     rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,...
                                                     identified_set_bounds,...
                                                     0.95);
data_output_folder = strcat( data_output_dir, dirname, 'Interacted_Moments/');
mkdir(data_output_folder);
save( strcat(data_output_folder, filename_graph,'_', 'rows_for_power_table'),...
    'lb_row_05', 'lb_row_50', 'lb_row_95',...
    'ub_row_05', 'ub_row_50', 'ub_row_95') ;                                              

%% Create excess length table
identified_set_length = identified_set_bounds(2) - identified_set_bounds(1);
identified_set_length_zerocutoff = identified_set_bounds_zerocutoff(2) - identified_set_bounds_zerocutoff(1);

%For the LF approaches, where we have a CS already, just take its length
%and subtract the identified set length
excess_lengths_lf = (confidence_sets_using_c_alpha(:,2) - confidence_sets_using_c_alpha(:,1)) - identified_set_length;
excess_lengths_lfn = (confidence_sets_using_c_lp_alpha(:,2) - confidence_sets_using_c_lp_alpha(:,1)) - identified_set_length;

%For the conditional and hybrid approaches, we compute the area by
%assigning the rejection value at each number to the closest point in the
%grid. The weight for each point is thus half the distance of the gap
%between the point on the left plus half the distance to the point on the
%right (the end points only get weight towards the interior, but this
%shouldn't matter since if set properly, the rejection probability is 1 at
%the endpoint)

gridWeightsVec = 0.5 * ([0,diff(beta0_grid)] + [diff(beta0_grid),0] );
gridWeightsMat = repmat( gridWeightsVec, size(full_rejection_grid_conditional,1), 1);

excess_lengths_conditional =  sum( gridWeightsMat .* (1 - full_rejection_grid_conditional ), 2) - identified_set_length;
excess_lengths_hybrid =  sum( gridWeightsMat .* (1 - full_rejection_grid_hybrid ), 2) - identified_set_length;

row_mean = cellfun( @(x) mean(x), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_05 = cellfun( @(x) quantile(x,.05), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_50 = cellfun( @(x) quantile(x,.50), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_95 = cellfun( @(x) quantile(x,.95), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});

save( strcat(data_output_folder, filename_graph,'_', 'rows_for_excess_length_table'),...
    'row_05', 'row_50', 'row_95', 'row_mean') ;                                              

row_mean_zerocutoff = row_mean - (identified_set_length_zerocutoff - identified_set_length);
row_05_zerocutoff = row_05 - (identified_set_length_zerocutoff - identified_set_length);
row_50_zerocutoff = row_50 - (identified_set_length_zerocutoff - identified_set_length);
row_95_zerocutoff = row_95 - (identified_set_length_zerocutoff - identified_set_length);

save( strcat(data_output_folder, filename_graph,'_', 'rows_for_excess_length_table_zerocutoff'),...
    'row_05_zerocutoff', 'row_50_zerocutoff', 'row_95_zerocutoff', 'row_mean_zerocutoff') ;                                              
%% Extract size at the upper and lower bound
upper_bound_index = find(beta0_grid == identified_set_bounds(2) );
lower_bound_index = find(beta0_grid == identified_set_bounds(1) );

id_set_indices = find(identified_set_bounds(1) <= beta0_grid & beta0_grid <= identified_set_bounds(2) );

upper_bound_sizes = cellfun( @(x) x(upper_bound_index), {rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,});
lower_bound_sizes = cellfun( @(x) x(lower_bound_index), {rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid});
                                                 
max_size_in_id_set = cellfun( @(x) max( x(id_set_indices)), {rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid});
                                                 
save( strcat(data_output_folder, filename_graph,'_', 'upper_and_lower_bound_size'),...
    'upper_bound_sizes', 'lower_bound_sizes', 'max_size_in_id_set') ;                                              

%% Extract size at the upper and lower bound using zerocutoff
upper_bound_index_zerocutoff = find(beta0_grid == identified_set_bounds_zerocutoff(2) );
lower_bound_index_zerocutoff = find(beta0_grid == identified_set_bounds_zerocutoff(1) );
id_set_indices_zerocutoff = find(identified_set_bounds_zerocutoff(1) <= beta0_grid & beta0_grid <= identified_set_bounds_zerocutoff(2) );


upper_bound_sizes_zerocutoff = cellfun( @(x) x(upper_bound_index_zerocutoff), {rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,});
lower_bound_sizes_zerocutoff = cellfun( @(x) x(lower_bound_index_zerocutoff), {rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid});
                                                 
max_size_in_id_set_zerocutoff = cellfun( @(x) max( x(id_set_indices_zerocutoff)), {rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid});
                                                                                                 
save( strcat(data_output_folder, filename_graph,'_', 'upper_and_lower_bound_size_zerocutoff'),...
    'upper_bound_sizes_zerocutoff', 'lower_bound_sizes_zerocutoff', 'max_size_in_id_set_zerocutoff') ;                                              


%%
display('Script complete');
