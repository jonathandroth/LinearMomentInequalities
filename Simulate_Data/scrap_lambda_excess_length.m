%% Create excess length table

full_rejection_grid_lf = lf_test_original;
full_rejection_grid_lfn = lf_test_modified;
full_rejection_grid_conditional = conditional_test;
full_rejection_grid_test = hybrid_test;


%For the conditional and hybrid approaches, we compute the area by
%assigning the rejection value at each number to the closest point in the
%grid. The weight for each point is thus half the distance of the gap
%between the point on the left plus half the distance to the point on the
%right (the end points only get weight towards the interior, but this
%shouldn't matter since if set properly, the rejection probability is 1 at
%the endpoint)

gridWeightsVec = 0.5 * ([0,diff(lambda_grid)] + [diff(lambda_grid),0] );
gridWeightsMat = repmat( gridWeightsVec, size(full_rejection_grid_conditional,1), 1);


acceptedAreaFn = @(grid) sum( gridWeightsMat .* (1 - grid ), 2);

identified_set_length = acceptedAreaFn(lambda_identified_set);

excessLengthFn = @(grid) acceptedAreaFn(grid) - identified_set_length ;

excess_lengths_lf = excessLengthFn(full_rejection_grid_lf);
excess_lengths_lfn = excessLengthFn(full_rejection_grid_lfn);
excess_lengths_conditional =  excessLengthFn(full_rejection_grid_conditional);
excess_lengths_hybrid =  excessLengthFn(full_rejection_grid_hybrid);



row_mean = cellfun( @(x) mean(x), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_05 = cellfun( @(x) quantile(x,.05), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_50 = cellfun( @(x) quantile(x,.50), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});
row_95 = cellfun( @(x) quantile(x,.95), {excess_lengths_lf, excess_lengths_lfn, excess_lengths_conditional, excess_lengths_hybrid});

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
    'row_05_with_infs', 'row_50_with_infs', 'row_95_with_infs', 'row_mean_with_infs') ;                                              

