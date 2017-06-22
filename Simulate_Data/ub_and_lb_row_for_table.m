%This function creates a row for a table with the lowerbounds for the
%identified set that has rejection probability equal to the argument power
%for each of the 4 methods, returned in lb_row.

%Likewise, ub_row returns the analgous upper_bound.

%Lowerbounds are calculated as being below the middle of the identified
%set. Upperbounds are calculated as being above the middle.

%Both rows contain the relevant bound for the identified set at the
%beginning as well

%The order is ID set boudh, LF, LFN, Cond, Hybrid

function [lb_row, ub_row] = ub_and_lb_row_for_table( beta0_grid, l_theta_grid,...
                                                     rejection_grid_c_alpha,...
                                                     rejection_grid_c_lp_alpha,...
                                                     rejection_grid_conditional,...
                                                     rejection_grid_hybrid,...
                                                     identified_set_bounds,...
                                                     power)

[lb_hybrid,ub_hybrid] = upper_and_lower_bounds_to_achieve_power( beta0_grid, rejection_grid_hybrid, identified_set_bounds,power);
[lb_conditional,ub_conditional] = upper_and_lower_bounds_to_achieve_power( beta0_grid, rejection_grid_conditional, identified_set_bounds,power);
[lb_lf,ub_lf] = upper_and_lower_bounds_to_achieve_power( l_theta_grid, rejection_grid_c_alpha, identified_set_bounds,power);
[lb_lfn,ub_lfn] = upper_and_lower_bounds_to_achieve_power( l_theta_grid, rejection_grid_c_lp_alpha, identified_set_bounds,power);

lb_row = [identified_set_bounds(1), lb_lf, lb_lfn, lb_conditional, lb_hybrid];
ub_row = [identified_set_bounds(2), ub_lf, ub_lfn, ub_conditional, ub_hybrid];


end