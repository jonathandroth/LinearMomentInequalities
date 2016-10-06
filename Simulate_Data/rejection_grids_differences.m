data_output_dir1 = '../../Output/Conditional_FullMatrix/Data/Calibrated_SigmaZeta/';
data_output_dirs2 = {'../../Output/Conditional_DiagonalMatrix/Data/Calibrated_SigmaZeta/',...
                     '../../Output/Conditional_OracleMatrix/Data/Calibrated_SigmaZeta/',...
                     '../../Output/UnconditionalMatrix/Data/Calibrated_SigmaZeta/'};
                 
figure_output_dirs = {'../../Figures/Differences/FullConditionalvDiagonalConditional/',...
                      '../../Figures/Differences/FullCondtionalvOracleConditional/',...
                      '../../Figures/Differences/FullConditionalvUnconditional/'};
                  
theta_g_grid = -150:5:100;
theta_c_grid = -250:10:510;
[ xgrid, ygrid] = meshgrid( theta_c_grid , theta_g_grid);
                  
                  
                  
for i = 1:length(data_output_dirs2)
    for moment_type = {'Basic_Moments/', 'Interacted_Moments/'}
    data_output_dir2 = data_output_dirs2(i);
    figure_output_dir = figure_output_dirs(i);

%data_output_dir2 = '../../Output/UnconditionalMatrix/Data/Calibrated_SigmaZeta/';
figure_output_dir = char(strcat(figure_output_dir, moment_type));
mkdir( figure_output_dir);

%Load  rejection grids cell for the two sets
 
rejection_grids_cell1 = load( char(strcat( data_output_dir1, moment_type, 'grid_cell.mat')) );

rejection_grids_cell2 = load( char(strcat( data_output_dir2, moment_type, 'grid_cell.mat')) );

if(strcmp(moment_type, 'Basic_Moments/'))
    rejection_grids_cell1 = rejection_grids_cell1.rejection_grids_cell;
    rejection_grids_cell2 = rejection_grids_cell2.rejection_grids_cell;
else
    rejection_grids_cell1 = rejection_grids_cell1.interacted_rejection_grids_cell;
    rejection_grids_cell2 = rejection_grids_cell2.interacted_rejection_grids_cell;
end



[rejection_prob_lf1, rejection_prob_rsw1, rejection_prob_conditional1, rejection_prob_hybrid1] = ...
    retrieve_rejection_grids_fn( rejection_grids_cell1, theta_c_grid, theta_g_grid );

[rejection_prob_lf2, rejection_prob_rsw2, rejection_prob_conditional2, rejection_prob_hybrid2] = ...
    retrieve_rejection_grids_fn( rejection_grids_cell2, theta_c_grid, theta_g_grid );


dif_lf = rejection_prob_lf1 - rejection_prob_lf2;
dif_rsw = rejection_prob_rsw1 - rejection_prob_rsw2;
dif_conditional = rejection_prob_conditional1 - rejection_prob_conditional2;
dif_hybrid = rejection_prob_hybrid1 - rejection_prob_hybrid2;

contour( xgrid, ygrid, dif_lf', [ [0,0.01,0.05,0.1,0.2], -[0.01,0.05,0.1,0.2] ], 'ShowText', 'on');
title('Difference in Rejection Probabilities - Least Favorable');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_output_dir, 'dif_rejection_probs_lf'), 'epsc');

contour( xgrid, ygrid, dif_rsw', [ [0,0.01,0.05,0.1,0.2], -[0.01,0.05,0.1,0.2] ], 'ShowText', 'on');
title('Difference in Rejection Probabilities - RSW');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_output_dir, 'dif_rejection_probs_rsw'), 'epsc');

contour( xgrid, ygrid, dif_conditional', [ [0,0.01,0.05,0.1,0.2], -[0.01,0.05,0.1,0.2] ], 'ShowText', 'on');
title('Difference in Rejection Probabilities - Conditional');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_output_dir, 'dif_rejection_probs_conditional'), 'epsc');


contour( xgrid, ygrid, dif_hybrid', [ [0,0.01,0.05,0.1,0.2], -[0.01,0.05,0.1,0.2] ], 'ShowText', 'on');
title('Difference in Rejection Probabilities - Hybrid');
xlabel('thetac');
ylabel('thetag');
saveas( gcf, strcat(figure_output_dir, 'dif_rejection_probs_hybrid'), 'epsc');

 end
end

