working_dir = '/Volumes/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data/';
cd(working_dir)

basic_inequalities_set_parameters
[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);

conditional_output_dir = '../../Output/Conditional_FullMatrix/Data/Calibrated_SigmaZeta/';
unconditional_output_dir = '../../Output/UnconditionalMatrix/Data/Calibrated_SigmaZeta/';

moment_type = 'Interacted_Moments/';
%moment_type = 'Basic_Moments/';

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
isoquant_value = 0.05;


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

saveas_tight( gcf, '../../Output/compare_rejection_probabilities_conditional', 'epsc');


%% 
%% Compare rejection probabilities for the 4 unconditional methods

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

saveas_tight( gcf, '../../Output/compare_rejection_probabilities_unconditional', 'epsc');

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

saveas_tight( gcf, '../../Output/compare_conditional_v_unconditional_lf', 'epsc');

contour( xgrid, ygrid, conditional_rejection_prob_rsw', [0,isoquant_value], 'ShowText','on', 'Color', colors(1,:) );
hold on;
contour( xgrid, ygrid, unconditional_rejection_prob_rsw', [0,isoquant_value], 'ShowText','on', 'Color', colors(2,:) );

%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN, 'Color',colors(1,:));
h(2) = plot(NaN,NaN, 'Color' ,colors(2,:));
legend(h, 'RSW (Conditional)', 'RSW (Unconditional');
hold off;

saveas_tight( gcf, '../../Output/compare_conditional_v_unconditional_rsw', 'epsc');



contour( xgrid, ygrid, conditional_rejection_prob_hybrid', [0,isoquant_value], 'ShowText','on', 'Color', colors(1,:) );
hold on;
contour( xgrid, ygrid, unconditional_rejection_prob_hybrid', [0,isoquant_value], 'ShowText','on', 'Color', colors(2,:) );

%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN, 'Color',colors(1,:));
h(2) = plot(NaN,NaN, 'Color' ,colors(2,:));
legend(h, 'Hybrid (Conditional)', 'Hybrid (Unconditional');
hold off;


saveas_tight( gcf, '../../Output/compare_conditional_v_unconditional_hybrid', 'epsc');



contour( xgrid, ygrid, conditional_rejection_prob_conditional', [0,isoquant_value], 'ShowText','on', 'Color', colors(1,:) );
hold on;
contour( xgrid, ygrid, unconditional_rejection_prob_conditional', [0,isoquant_value], 'ShowText','on', 'Color', colors(2,:) );

%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN, 'Color',colors(1,:));
h(2) = plot(NaN,NaN, 'Color' ,colors(2,:));
legend(h, 'Conditional (Conditional Var)', 'Conditional (Unconditional Var)');
hold off;


saveas_tight( gcf, '../../Output/compare_conditional_v_unconditional_conditional', 'epsc');

%% Create a table with average volume of accepted region
dim1 = theta_g_grid(2) - theta_g_grid(1);
dim2 = theta_c_grid(2) - theta_c_grid(1);

conditional_areas = cellfun( @(grid) calculate_accepted_area(grid,dim1,dim2) , {conditional_rejection_prob_lf;...
                                     conditional_rejection_prob_rsw;...
                                     conditional_rejection_prob_hybrid;...
                                     conditional_rejection_prob_conditional} )

unconditional_areas = cellfun( @(grid) calculate_accepted_area(grid,dim1,dim2) , {unconditional_rejection_prob_lf;...
                                     unconditional_rejection_prob_rsw;...
                                     unconditional_rejection_prob_hybrid;...
                                     unconditional_rejection_prob_conditional} )
                                 
fh = fopen('../../Output/compare_grid_areas.tex','w');                                 
diary '../../Output/compare_grid_areas.tex'
 latex2( '%0.3g', {'', 'Conditional Var.', 'Unconditional Var.'},...
     {'LF';'RSW';'Hybrid'; 'Conditional'} ,conditional_areas, unconditional_areas );
 diary off
