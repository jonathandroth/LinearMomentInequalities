working_dir = '/Volumes/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data';
cd(working_dir)

basic_inequalities_set_parameters
[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);

conditional_output_dir = '../../Output/Conditional_FullMatrix/Data/Calibrated_SigmaZeta/';
unconditional_output_dir = '../../Output/UnconditionalMatrix/Data/Calibrated_SigmaZeta/';

%moment_type = 'Interacted_Moments/';
moment_type = 'Basic_Moments/';

conditional_rejection_grids_cell = load( char(strcat( conditional_output_dir, moment_type, 'grid_cell.mat')) );

unconditional_rejection_grids_cell = load( char(strcat( unconditional_output_dir, moment_type, 'grid_cell.mat')) );


if(strcmp(moment_type, 'Basic_Moments/'))
    conditional_rejection_grids_cell = conditional_rejection_grids_cell.rejection_grids_cell;
    unconditional_rejection_grids_cell = unconditional_rejection_grids_cell.rejection_grids_cell;
else
    conditional_rejection_grids_cell = conditional_rejection_grids_cell.interacted_rejection_grids_cell;
    unconditional_rejection_grids_cell = unconditional_rejection_grids_cell.interacted_rejection_grids_cell;
end



[conditional_rejection_prob_lf, conditional_rejection_prob_rsw, conditional_rejection_prob_conditional, conditional_rejection_prob_hybrid] = ...
    retrieve_rejection_grids_fn( conditional_rejection_grids_cell, theta_c_grid, theta_g_grid );

[unconditional_rejection_prob_lf, unconditional_rejection_prob_rsw, unconditional_rejection_prob_conditional, unconditional_rejection_prob_hybrid] = ...
    retrieve_rejection_grids_fn( unconditional_rejection_grids_cell, theta_c_grid, theta_g_grid );


%% Compare rejection probabilities for the 4 conditional methods

%We will show one line for each method. This specifies the isoquant value
%for the rejection prob
isoquant_value = 0.5;


colors = get(gca,'colororder');

contour( xgrid, ygrid, conditional_rejection_prob_lf', [0,isoquant_value], 'ShowText','on', 'Color', colors(1,:) );
hold on;
contour( xgrid, ygrid, conditional_rejection_prob_rsw', [0,isoquant_value],'ShowText','on', 'Color', colors(2,:) );
contour( xgrid, ygrid, conditional_rejection_prob_hybrid', [0,isoquant_value], 'ShowText','on', 'Color', colors(3,:));
%contour( xgrid, ygrid, conditional_rejection_prob_conditional', [0,isoquant_value], 'ShowText','on','Color', 'g');

h = zeros(3, 1);
h(1) = plot(NaN,NaN, 'Color',colors(1,:));
h(2) = plot(NaN,NaN, 'Color' ,colors(2,:));
h(3) = plot(NaN,NaN, 'Color' ,colors(3,:));
legend(h, 'LF', 'RSW', 'Hybrid');

hold off;


%% 
%% Compare rejection probabilities for the 4 unconditional methods

%We will show one line for each method. This specifies the isoquant value
%for the rejection prob
isoquant_value = 0.1;


colors = get(gca,'colororder');

contour( xgrid, ygrid, unconditional_rejection_prob_lf', [0,isoquant_value], 'ShowText','on', 'Color', colors(1,:) );
hold on;
contour( xgrid, ygrid, unconditional_rejection_prob_rsw', [0,isoquant_value],'ShowText','on', 'Color', colors(2,:) );
contour( xgrid, ygrid, unconditional_rejection_prob_hybrid', [0,isoquant_value], 'ShowText','on', 'Color', colors(3,:));
%contour( xgrid, ygrid, unconditional_rejection_prob_conditional', [0,isoquant_value], 'ShowText','on','Color', 'g');

h = zeros(3, 1);
h(1) = plot(NaN,NaN, 'Color',colors(1,:));
h(2) = plot(NaN,NaN, 'Color' ,colors(2,:));
h(3) = plot(NaN,NaN, 'Color' ,colors(3,:));
legend(h, 'LF', 'RSW', 'Hybrid');

hold off;


%% Compare conditional versus unconditional

contour( xgrid, ygrid, conditional_rejection_prob_lf', [0,isoquant_value], 'ShowText','on', 'Color', colors(1,:) );
hold on;
contour( xgrid, ygrid, unconditional_rejection_prob_lf', [0,isoquant_value], 'ShowText','on', 'Color', colors(2,:) );

%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN, 'Color',colors(1,:));
h(2) = plot(NaN,NaN, 'Color' ,colors(2,:));
legend(h, 'LF (Conditional)', 'LF (Unconditional');
hold off;


contour( xgrid, ygrid, conditional_rejection_prob_rsw', [0,isoquant_value], 'ShowText','on', 'Color', colors(1,:) );
hold on;
contour( xgrid, ygrid, unconditional_rejection_prob_rsw', [0,isoquant_value], 'ShowText','on', 'Color', colors(2,:) );

%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN, 'Color',colors(1,:));
h(2) = plot(NaN,NaN, 'Color' ,colors(2,:));
legend(h, 'RSW (Conditional)', 'RSW (Unconditional');
hold off;



contour( xgrid, ygrid, conditional_rejection_prob_hybrid', [0,isoquant_value], 'ShowText','on', 'Color', colors(1,:) );
hold on;
contour( xgrid, ygrid, unconditional_rejection_prob_hybrid', [0,isoquant_value], 'ShowText','on', 'Color', colors(2,:) );

%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN, 'Color',colors(1,:));
h(2) = plot(NaN,NaN, 'Color' ,colors(2,:));
legend(h, 'Hybrid (Conditional)', 'Hybrid (Unconditional');
hold off;