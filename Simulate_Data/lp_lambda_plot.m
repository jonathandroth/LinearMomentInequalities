 
 data_output_folder = strcat( data_output_dir, dirname, 'Interacted_Moments/');

 ds_dir = strcat( data_output_dir, 'Interacted_Moments/');
 load(  strcat(ds_dir, 'lambda_identified_set') );
 load(  strcat(ds_dir, 'lambda_results') );

addpath('./export-fig')

getrow = @(x,row) x(row,:);

clf
hold('on');

pbaspect([1.2 1 1])
l1 = plot(lambda_vec , mean(hybrid_test), 'Color', getrow( get(gca,'colororder'),1 ), 'LineStyle', '-' )
l2 = plot(lambda_vec , mean(conditional_test), 'Color', getrow( get(gca,'colororder'),5 ), 'LineStyle', ':')
l3 = plot(lambda_vec, mean(lf_test_modified), 'Color', getrow( get(gca,'colororder'),3 ), 'LineStyle', '-.')

 
identified_set_max = max( lambda_vec( lambda_identified_set == 1) );
identified_set_min = min( lambda_vec( lambda_identified_set == 1) );


if(~isempty(identified_set_max) )
    l5 = plot( [identified_set_max; identified_set_max], [0;1], 'LineStyle', '--', 'Color',  'r');
else
    warning('Didnt find any lambdas in identfied set');
end

if(~isempty(identified_set_min) && identified_set_min ~= min(lambda_vec)  )
   l6 = plot( [identified_set_min; identified_set_min], [0;1], 'LineStyle', '--', 'Color',  'r');
end

legend( 'Hybrid', 'Conditional', 'LF', 'Identified Set Bound', 'Location',....
    'southoutside', 'Orientation', 'horizontal' );
ylabel('Rejection Probability');
%xlabel('Beta');
%title('Rejection Probabilities for Beta');

set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf, 'Type', 'Line'),'LineWidth',4); %Linewidth for plot lines


%saveas( gcf, strcat(figures_output_dir,filename_graph ), 'epsc');
export_fig(strcat(figures_output_dir,filename_graph,'.pdf'));

clf
display( strcat(figures_output_dir,filename_graph ) )

%% Make plot for sCC and sRCC test and LFP


hold('on');

pbaspect([1.4 1 1])
l1 = plot(lambda_vec , mean(hybrid_test), 'Color', getrow( get(gca,'colororder'),1 ), 'LineStyle', '-' )
l2 = plot(lambda_vec , mean(lf_test_original), 'Color', getrow( get(gca,'colororder'),3 ), 'LineStyle', '-.' )

l3 = plot(lambda_vec , mean(cc_test), 'Color', getrow( get(gca,'colororder'),4 ), 'LineStyle', '-.' )
l4 = plot(lambda_vec , mean(rcc_test), 'Color', getrow( get(gca,'colororder'),5 ), 'LineStyle', ':' )
 
identified_set_max = max( lambda_vec( lambda_identified_set == 1) );
identified_set_min = min( lambda_vec( lambda_identified_set == 1) );


if(~isempty(identified_set_max) )
    l5 = plot( [identified_set_max; identified_set_max], [0;1], 'LineStyle', '--', 'Color',  'r');
else
    warning('Didnt find any lambdas in identfied set');
end

if(~isempty(identified_set_min) && identified_set_min ~= min(lambda_vec)  )
   l6 = plot( [identified_set_min; identified_set_min], [0;1], 'LineStyle', '--', 'Color',  'r');
end

legend( 'Hybrid', 'LFP', 'sCC', 'sRCC', 'ID Set Bound', 'Location',....
    'southoutside', 'Orientation', 'horizontal' );
ylabel('Rejection Probability');
%xlabel('Beta');
%title('Rejection Probabilities for Beta');

set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf, 'Type', 'Line'),'LineWidth',4); %Linewidth for plot lines


%saveas( gcf, strcat(figures_output_dir,filename_graph ), 'epsc');
export_fig(strcat(figures_output_dir,filename_graph,'_compare_to_Cox_and_Shi','.pdf'));

clf
%% Make rows for excess length table

full_rejection_grid_lf = lf_test_original;
full_rejection_grid_lfn = lf_test_modified;
full_rejection_grid_conditional = conditional_test;
full_rejection_grid_hybrid = hybrid_test;

full_rejection_grid_cc = cc_test;
full_rejection_grid_rcc = rcc_test;

rejection_grid_lf = mean(lf_test_original,1);
rejection_grid_lfn = mean(lf_test_modified,1);
rejection_grid_conditional = mean(conditional_test,1);
rejection_grid_hybrid = mean(hybrid_test,1);

rejection_grid_cc = mean(cc_test,1);
rejection_grid_rcc = mean(rcc_test,1);

%For the conditional and hybrid approaches, we compute the area by
%assigning the rejection value at each number to the closest point in the
%grid. The weight for each point is thus half the distance of the gap
%between the point on the left plus half the distance to the point on the
%right (the end points only get weight towards the interior, but this
%shouldn't matter since if set properly, the rejection probability is 1 at
%the endpoint)

gridWeightsVec = 0.5 * ([0,diff(lambda_vec')] + [diff(lambda_vec'),0] );
gridWeightsMat = repmat( gridWeightsVec, size(full_rejection_grid_conditional,1), 1);



identified_set_length = sum( gridWeightsVec .* lambda_identified_set');
identified_set_length_zerocutoff = sum( gridWeightsVec .* lambda_identified_set_zerocutoff');

acceptedAreaFn = @(grid) sum( gridWeightsMat .* (1 - grid ), 2);
excessLengthFn = @(grid) acceptedAreaFn(grid) - identified_set_length ;

excess_lengths_lf = excessLengthFn(full_rejection_grid_lf);
excess_lengths_lfn = excessLengthFn(full_rejection_grid_lfn);
excess_lengths_conditional =  excessLengthFn(full_rejection_grid_conditional);
excess_lengths_hybrid =  excessLengthFn(full_rejection_grid_hybrid);

excess_lengths_cc =  sum( gridWeightsMat .* (1 - full_rejection_grid_cc ), 2) - identified_set_length;
excess_lengths_rcc =  sum( gridWeightsMat .* (1 - full_rejection_grid_rcc ), 2) - identified_set_length;


row_mean = cellfun( @(x) mean(x), {excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid, excess_lengths_lf});
row_05 = cellfun( @(x) quantile(x,.05), {excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid, excess_lengths_lf});
row_50 = cellfun( @(x) quantile(x,.50), {excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid, excess_lengths_lf});
row_95 = cellfun( @(x) quantile(x,.95), {excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid, excess_lengths_lf});


row_mean_coxandshi = cellfun( @(x) mean(x), {excess_lengths_cc, excess_lengths_rcc});
row_05_coxandshi = cellfun( @(x) quantile(x,.05), {excess_lengths_cc, excess_lengths_rcc});
row_50_coxandshi = cellfun( @(x) quantile(x,.50), {excess_lengths_cc, excess_lengths_rcc});
row_95_coxandshi = cellfun( @(x) quantile(x,.95), {excess_lengths_cc, excess_lengths_rcc});


row_mean_zerocutoff = row_mean - (identified_set_length_zerocutoff - identified_set_length);
row_05_zerocutoff = row_05 - (identified_set_length_zerocutoff - identified_set_length);
row_50_zerocutoff = row_50 - (identified_set_length_zerocutoff - identified_set_length);
row_95_zerocutoff = row_95 - (identified_set_length_zerocutoff - identified_set_length);

row_mean_zerocutoff_coxandshi = row_mean_coxandshi - (identified_set_length_zerocutoff - identified_set_length);
row_05_zerocutoff_coxandshi = row_05_coxandshi - (identified_set_length_zerocutoff - identified_set_length);
row_50_zerocutoff_coxandshi = row_50_coxandshi - (identified_set_length_zerocutoff - identified_set_length);
row_95_zerocutoff_coxandshi = row_95_coxandshi - (identified_set_length_zerocutoff - identified_set_length);



%Compute the various statistics treating unbounded instances as infinity,
%rather than truncated at the endpoint, as done above
unbounded_lf = full_rejection_grid_lf(:,end) < 1;
unbounded_lfn = full_rejection_grid_lfn(:,end) < 1;
unbounded_conditional = full_rejection_grid_conditional(:,end) < 1;
unbounded_hybrid = full_rejection_grid_hybrid(:,end) < 1;

unbounded_cc = full_rejection_grid_cc(:,end) < 1;
unbounded_rcc = full_rejection_grid_rcc(:,end) < 1;

excess_lengths_lf_with_infs = ifelse( unbounded_lf, Inf , excess_lengths_lf);
excess_lengths_lfn_with_infs = ifelse( unbounded_lfn, Inf , excess_lengths_lfn);
excess_lengths_conditional_with_infs = ifelse( unbounded_conditional, Inf , excess_lengths_conditional);
excess_lengths_hybrid_with_infs = ifelse( unbounded_hybrid, Inf , excess_lengths_hybrid);

excess_lengths_cc_with_infs = ifelse( unbounded_cc, Inf , excess_lengths_cc);
excess_lengths_rcc_with_infs = ifelse( unbounded_rcc, Inf , excess_lengths_rcc);


row_mean_with_infs = cellfun( @(x) mean(x), {excess_lengths_lfn_with_infs, excess_lengths_conditional_with_infs, excess_lengths_hybrid_with_infs, excess_lengths_lf_with_infs});
row_05_with_infs = cellfun( @(x) quantile(x,.05), {excess_lengths_lfn_with_infs, excess_lengths_conditional_with_infs, excess_lengths_hybrid_with_infs, excess_lengths_lf_with_infs});
row_50_with_infs = cellfun( @(x) quantile(x,.50), {excess_lengths_lfn_with_infs, excess_lengths_conditional_with_infs, excess_lengths_hybrid_with_infs, excess_lengths_lf_with_infs});
row_95_with_infs = cellfun( @(x) quantile(x,.95), {excess_lengths_lfn_with_infs, excess_lengths_conditional_with_infs, excess_lengths_hybrid_with_infs, excess_lengths_lf_with_infs});


row_mean_with_infs_coxandshi = cellfun( @(x) mean(x), {excess_lengths_cc_with_infs, excess_lengths_rcc_with_infs});
row_05_with_infs_coxandshi = cellfun( @(x) quantile(x,.05), {excess_lengths_cc_with_infs, excess_lengths_rcc_with_infs});
row_50_with_infs_coxandshi = cellfun( @(x) quantile(x,.50), {excess_lengths_cc_with_infs, excess_lengths_rcc_with_infs});
row_95_with_infs_coxandshi = cellfun( @(x) quantile(x,.95), {excess_lengths_cc_with_infs, excess_lengths_rcc_with_infs});


save( strcat(data_output_folder, filename_graph,'_', 'rows_for_excess_length_table_lambda'),...
    'row_05', 'row_50', 'row_95', 'row_mean',...
    'row_05_zerocutoff', 'row_50_zerocutoff', 'row_95_zerocutoff', 'row_mean_zerocutoff',...
    'row_05_with_infs', 'row_50_with_infs', 'row_95_with_infs', 'row_mean_with_infs',...
    'row_05_coxandshi', 'row_50_coxandshi', 'row_95_coxandshi', 'row_mean_coxandshi',...
    'row_05_zerocutoff_coxandshi', 'row_50_zerocutoff_coxandshi', 'row_95_zerocutoff_coxandshi', 'row_mean_zerocutoff_coxandshi',...
    'row_05_with_infs_coxandshi', 'row_50_with_infs_coxandshi', 'row_95_with_infs_coxandshi', 'row_mean_with_infs_coxandshi') ;                                              



%% Make row for size table

upper_bound_index = max( find(lambda_identified_set == 1 ) );
id_set_indices = find(lambda_identified_set == 1 );

upper_bound_sizes = cellfun( @(x) x(upper_bound_index), {...
                                                     rejection_grid_lfn,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,...
                                                     rejection_grid_lf});

                                                 
max_size_in_id_set = cellfun( @(x) max( x(id_set_indices)), {...
                                                     rejection_grid_lfn,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,...
                                                     rejection_grid_lf});

max_size_in_id_set_coxandshi = cellfun( @(x) max( x(id_set_indices)), ...
                                                    {rejection_grid_cc,...
                                                     rejection_grid_rcc});     
                                                 
save( strcat(data_output_folder, filename_graph,'_', 'upper_bound_size_lambda'),...
    'upper_bound_sizes', 'max_size_in_id_set', 'max_size_in_id_set_coxandshi') ;                                              


%% Make row for size table - zerocutoff
upper_bound_index_zerocutoff = max( find(lambda_identified_set_zerocutoff == 1 ) );
id_set_indices_zerocutoff = find(lambda_identified_set_zerocutoff == 1 );


upper_bound_sizes_zerocutoff = cellfun( @(x) x(upper_bound_index_zerocutoff), {...
                                                     rejection_grid_lfn,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,...
                                                     rejection_grid_lf});

max_size_in_id_set_zerocutoff = cellfun( @(x) max( x(id_set_indices_zerocutoff)), {...
                                                     rejection_grid_lfn,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,...
                                                     rejection_grid_lf});
                                                 
max_size_in_id_set_zerocutoff_coxandshi = cellfun( @(x) max( x(id_set_indices_zerocutoff)), ...
                                                    {rejection_grid_cc,...
                                                     rejection_grid_rcc});                                                 
                                                 
                                                 
save( strcat(data_output_folder, filename_graph,'_', 'upper_bound_size_lambda_zerocutoff'),...
    'upper_bound_sizes_zerocutoff', 'max_size_in_id_set_zerocutoff','max_size_in_id_set_zerocutoff_coxandshi') ; 



%% Make row for timing table

if(exist('timing_vec_lf_lambda'))
timing_row_allcombined = ...
              [nanmean(sum(timing_vec_lf_lambda,2)),...
              nanmean(sum(timing_vec_conditional_lambda,2)),...
              nanmean(sum(timing_vec_hybrid_lambda,2)),...
              nanmean(sum(timing_vec_lfp_lambda,2)),...
              nanmean(sum(timing_vec_cc_lambda,2)),...
              nanmean(sum(timing_vec_rcc_lambda,2))...
              ];
save( strcat(data_output_folder, filename_graph,'_', 'timing_row_allcombined'),...
    'timing_row_allcombined') ;                                              
          
end
