 
 data_output_folder = strcat( data_output_dir, dirname, 'Interacted_Moments/');

 ds_dir = strcat( data_output_dir, 'Interacted_Moments/');
 load(  strcat(ds_dir, 'lambda_identified_set') );
 load(  strcat(ds_dir, 'lambda_results') );
 
  
 plot(lambda_vec, [mean(lf_test_original);...
                  mean(lf_test_modified);...
                  mean(conditional_test);...
                  mean(hybrid_test)]) 


identified_set_max = max( lambda_vec( lambda_identified_set == 1) );
identified_set_min = min( lambda_vec( lambda_identified_set == 1) );


if(~isempty(identified_set_max) )
    line( [identified_set_max; identified_set_max], [0;1], 'LineStyle', '--', 'Color',  'r');
else
    warning('Didnt find any lambdas in identfied set');
end

if(~isempty(identified_set_min) && identified_set_min ~= min(lambda_vec)  )
    line( [identified_set_min; identified_set_min], [0;1], 'LineStyle', '--', 'Color',  'r');
end

legend( 'LFP','LF', 'Conditional', 'Hybrid', 'Identified Set Bound', 'Location','eastoutside' );
ylabel('Rejection Probability');
xlabel('Beta');
%title('Rejection Probabilities for Beta');

saveas( gcf, strcat(figures_output_dir,filename_graph ), 'epsc');
%export_fig(strcat(figures_output_dir,filename_graph,'.pdf'));

display( strcat(figures_output_dir,filename_graph ) )


%% Make rows for excess length table

full_rejection_grid_lf = lf_test_original;
full_rejection_grid_lfn = lf_test_modified;
full_rejection_grid_conditional = conditional_test;
full_rejection_grid_hybrid = hybrid_test;

rejection_grid_lf = mean(lf_test_original,1);
rejection_grid_lfn = mean(lf_test_modified,1);
rejection_grid_conditional = mean(conditional_test,1);
rejection_grid_hybrid = mean(hybrid_test,1);

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



row_mean = cellfun( @(x) mean(x), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_05 = cellfun( @(x) quantile(x,.05), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_50 = cellfun( @(x) quantile(x,.50), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_95 = cellfun( @(x) quantile(x,.95), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});


row_mean_zerocutoff = row_mean - (identified_set_length_zerocutoff - identified_set_length);
row_05_zerocutoff = row_05 - (identified_set_length_zerocutoff - identified_set_length);
row_50_zerocutoff = row_50 - (identified_set_length_zerocutoff - identified_set_length);
row_95_zerocutoff = row_95 - (identified_set_length_zerocutoff - identified_set_length);

%Compute the various statistics treating unbounded instances as infinity,
%rather than truncated at the endpoint, as done above
unbounded_lf = full_rejection_grid_lf(:,end) < 1;
unbounded_lfn = full_rejection_grid_lfn(:,end) < 1;
unbounded_conditional = full_rejection_grid_conditional(:,end) < 1;
unbounded_hybrid = full_rejection_grid_hybrid(:,end) < 1;

excess_lengths_lf_with_infs = ifelse( unbounded_lf, Inf , excess_lengths_lf);
excess_lengths_lfn_with_infs = ifelse( unbounded_lfn, Inf , excess_lengths_lfn);
excess_lengths_conditional_with_infs = ifelse( unbounded_conditional, Inf , excess_lengths_conditional);
excess_lengths_hybrid_with_infs = ifelse( unbounded_hybrid, Inf , excess_lengths_hybrid);

row_mean_with_infs = cellfun( @(x) mean(x), {excess_lengths_lf_with_infs, excess_lengths_lfn_with_infs, excess_lengths_conditional_with_infs, excess_lengths_hybrid_with_infs});
row_05_with_infs = cellfun( @(x) quantile(x,.05), {excess_lengths_lf_with_infs, excess_lengths_lfn_with_infs, excess_lengths_conditional_with_infs, excess_lengths_hybrid_with_infs});
row_50_with_infs = cellfun( @(x) quantile(x,.50), {excess_lengths_lf_with_infs, excess_lengths_lfn_with_infs, excess_lengths_conditional_with_infs, excess_lengths_hybrid_with_infs});
row_95_with_infs = cellfun( @(x) quantile(x,.95), {excess_lengths_lf_with_infs, excess_lengths_lfn_with_infs, excess_lengths_conditional_with_infs, excess_lengths_hybrid_with_infs});



save( strcat(data_output_folder, filename_graph,'_', 'rows_for_excess_length_table_lambda'),...
    'row_05', 'row_50', 'row_95', 'row_mean',...
    'row_05_zerocutoff', 'row_50_zerocutoff', 'row_95_zerocutoff', 'row_mean_zerocutoff',...
    'row_05_with_infs', 'row_50_with_infs', 'row_95_with_infs', 'row_mean_with_infs') ;                                              


%% Make row for size table

upper_bound_index = max( find(lambda_identified_set == 1 ) );
id_set_indices = find(lambda_identified_set == 1 );

upper_bound_sizes = cellfun( @(x) x(upper_bound_index), {rejection_grid_lf,...
                                                     rejection_grid_lfn,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,});

                                                 
max_size_in_id_set = cellfun( @(x) max( x(id_set_indices)), {rejection_grid_lf,...
                                                     rejection_grid_lfn,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid});

                                                 
save( strcat(data_output_folder, filename_graph,'_', 'upper_bound_size_lambda'),...
    'upper_bound_sizes', 'max_size_in_id_set') ;                                              


%% Make row for size table - zerocutoff
upper_bound_index_zerocutoff = max( find(lambda_identified_set_zerocutoff == 1 ) );
id_set_indices_zerocutoff = find(lambda_identified_set_zerocutoff == 1 );


upper_bound_sizes_zerocutoff = cellfun( @(x) x(upper_bound_index_zerocutoff), {rejection_grid_lf,...
                                                     rejection_grid_lfn,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,});

max_size_in_id_set_zerocutoff = cellfun( @(x) max( x(id_set_indices_zerocutoff)), {rejection_grid_lf,...
                                                     rejection_grid_lfn,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid});
                                                 
save( strcat(data_output_folder, filename_graph,'_', 'upper_bound_size_lambda_zerocutoff'),...
    'upper_bound_sizes_zerocutoff', 'max_size_in_id_set_zerocutoff') ;                                              
