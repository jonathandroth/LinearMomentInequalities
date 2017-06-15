clf
addpath('./BreakXAxis_2')

plot( repmat( l_theta_grid', 1, 2), [rejection_grid_c_alpha, rejection_grid_c_lp_alpha ]);

line( beta0_grid , rejection_grid_conditional, 'Color', getrow( get(gca,'colororder'),3 ) )
line( beta0_grid , rejection_grid_hybrid, 'Color', getrow( get(gca,'colororder'),4 ) )

line( [identified_set_bounds(1);identified_set_bounds(1)] ,[0;1], 'LineStyle', '--', 'Color',  'r');
line( [identified_set_bounds(2);identified_set_bounds(2)] ,[0;1], 'LineStyle', '--', 'Color',  'r');


legend( 'LF','LFN', 'Conditional', 'Hybrid', 'Identified Set Boundary',  'Location','eastoutside' );
ylabel('Rejection Probability');


if( exist('xlim_graph') ==1)
    %If manual bounds are specified for the x-axis limit, impose these
    xlim( xlim_graph )    
    
    %Impose tick width if specified; otherwise 5
    if(~exist('xtick_width'))
        xtick_width = 5;
    end
    set(gca,'XTick',[xlim_graph(1):xtick_width:xlim_graph(2)])
    
    %Create a break in the x-axis
    if(exist('xsplit_graph'))
        %breakxaxis(xsplit_graph)
        BreakXAxis(beta0_grid, rejection_grid_conditional, xsplit_graph(1),xsplit_graph(2), 1)
    end


end

    





