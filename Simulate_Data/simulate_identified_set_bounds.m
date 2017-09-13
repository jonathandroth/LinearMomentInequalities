tic;
%1thetac_interacted params

%specify whether to combine the moments for theta_g, or to have separate
%ones for each firm
combine_theta_g_moments = 1;


%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/1Thetac/Interacted_Moments/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/LP_figures/1Thetac/Interacted_Moments/';


%Specify the groups of firms that have different coefficients
F_group_cell_moments = {[1;2;3;4;5;6;7;8;9]};


%num_F_groups_moments = size(F_group_cell_moments,1);
%l = [1; zeros(num_F_groups_moments-1,1); mean_g];

%Run the main script (this does confidence sets for the mean weight)
if( exist('xlim_graph_meanweight'))
    xlim_graph = xlim_graph_meanweight;
else
    xlim_graph = [-150;175];
end

if( exist('xsplit_graph_meanweight'))
   xsplit_graph = xsplit_graph_meanweight;
end

filename_graph = 'Mean_Weight_Rejection_Probabilities';

lp_set_parameters
%%
load('../../Output/Simulated_Data/all_simulation_params')


numSimulations = 100;

boundsMat_meanweight = NaN(numSimulations,2);
ub_slackCell_meanweight = cell(numSimulations,1);
lb_slackCell_meanweight = cell(numSimulations,1);
boundsMat_thetag = NaN(numSimulations,2);
ub_slackCell_thetag = cell(numSimulations,1);
lb_slackCell_thetag = cell(numSimulations,1);


for(sim = 1:numSimulations)

sim

[J_t_array, J_tminus1_array, Pi_array, F_array, G_array, Pi_star_array, Eta_jt_shocks_array, Eta_t_vec] = simulate_data(sim, F , J ,T, burnout, sigma_nu, sigma_eps, sigma_w, sigma_zetaj, sigma_zetajft, rho, lambda,theta_c,theta_g, g_vec, mu_f);

[y_T, X_T] = lp_create_moments_for_identified_set_fn(F_group_cell_moments, F_group_cell_parameters, F_array, G_array, Eta_jt_shocks_array, Eta_t_vec, Pi_array, J_t_array, J_tminus1_array, use_basic_moments, lambda, combine_theta_g_moments);

%Compute ID set for meanweight
mean_g = mean(G_array(1,:,1));
l_meanweight = [ones(num_F_groups_parameters,1) / num_F_groups_parameters; mean_g];
[identified_set_bounds_for_meanweight, slack_lb_meanweight, slack_ub_meanweight] = cs_linear_delta_lp_fn(y_T,X_T,l_meanweight,0);


boundsMat_meanweight(sim,:) = identified_set_bounds_for_meanweight;
ub_slackCell_meanweight{sim} = slack_ub_meanweight;
lb_slackCell_meanweight{sim} = slack_lb_meanweight;


l_thetag = [ zeros( num_F_groups_parameters ,1) ; 1 ];
[identified_set_bounds_for_thetag, slack_lb_thetag, slack_ub_thetag] = cs_linear_delta_lp_fn(y_T,X_T,l_thetag,0);


boundsMat_thetag(sim,:) = identified_set_bounds_for_thetag;
ub_slackCell_thetag{sim} = slack_ub_thetag;
lb_slackCell_thetag{sim} = slack_lb_thetag;

end

%%
%Create functions to convert the ub and lb cells to a vector that says what
%fraction of the time the moments are binding at the optimum
convertToMat = @(c) cell2mat(c');
approxZero = @(x) (abs(x) <= 10^(-6) );
summariseBindingMoments = @(c) mean( approxZero( convertToMat(c) ) , 2); 

ub_bindingMoments_meanweight = summariseBindingMoments(ub_slackCell_meanweight);
lb_bindingMoments_meanweight = summariseBindingMoments(lb_slackCell_meanweight);

ub_bindingMoments_thetag = summariseBindingMoments(ub_slackCell_thetag);
lb_bindingMoments_thetag = summariseBindingMoments(lb_slackCell_thetag);


save( strcat(data_output_dir,'identified_set_results.mat'), ...
    'ub_bindingMoments_meanweight', 'ub_bindingMoments_thetag',...
    'lb_bindingMoments_meanweight', 'lb_bindingMoments_thetag',...
    'ub_slackCell_meanweight', 'ub_slackCell_thetag',...
    'lb_slackCell_meanweight', 'lb_slackCell_thetag',...
    'boundsMat_meanweight', 'boundsMat_thetag')


addpath('./swtest')
addpath('./export-fig')
addpath('./matrix2latex')
%%Histograms of bounds
%%Graphs of q-q plots with Shapiro-Wilk Test for Normality
clf
hist(boundsMat_thetag(:,1) );
title('Lower Bound for Thetag');
export_fig(strcat(figures_output_dir, 'histogram_lb_thetag.pdf') )
clf

hist(boundsMat_thetag(:,2) );
title('Upper Bound for Thetag');
export_fig(strcat(figures_output_dir, 'histogram_ub_thetag.pdf') )
clf

hist(boundsMat_meanweight(:,1) );
title('Lower Bound for Mean-Weight Cost');
export_fig(strcat(figures_output_dir, 'histogram_lb_meanweight.pdf') )
clf

hist(boundsMat_meanweight(:,2) );
title('Upper Bound for Mean-Weight Cost');
export_fig(strcat(figures_output_dir, 'histogram_ub_meanweight.pdf') )
clf


%%Graphs of q-q plots with Shapiro-Wilk Test for Normality
clf
qqplot_with_pvalue(boundsMat_thetag(:,1) );
title('Lower Bound for Thetag');
export_fig(strcat(figures_output_dir, 'qqplot_lb_thetag.pdf') )
clf

qqplot_with_pvalue(boundsMat_thetag(:,2) );
title('Upper Bound for Thetag');
export_fig(strcat(figures_output_dir, 'qqplot_ub_thetag.pdf') )
clf

qqplot_with_pvalue(boundsMat_meanweight(:,1) );
title('Lower Bound for Mean-Weight Cost');
export_fig(strcat(figures_output_dir, 'qqplot_lb_meanweight.pdf') )
clf

qqplot_with_pvalue(boundsMat_meanweight(:,2) );
title('Upper Bound for Mean-Weight Cost');
export_fig(strcat(figures_output_dir, 'qqplot_ub_meanweight.pdf') )
clf


matrix2latex( [(1:size(lb_bindingMoments_thetag,1))',...
               lb_bindingMoments_thetag, ub_bindingMoments_thetag , ...
               lb_bindingMoments_meanweight, ub_bindingMoments_meanweight] ,...
               strcat( figures_output_dir, "binding_moments.tex" ),...
               'columnLabels', {'Moment \#','LB Thetag', 'UB Thetag', 'LB Meanweight', 'UB Meanweight'});

toc;
