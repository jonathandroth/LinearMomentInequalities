function [lower_bound, upper_bound] = upper_and_lower_bounds_to_achieve_power( thetagrid, rejectiongrid, identified_set_bounds,power)

middle_of_identified_set = mean( identified_set_bounds);

upper_theta_indices = find( thetagrid >= middle_of_identified_set);
lower_theta_indices = find( thetagrid <= middle_of_identified_set);

lower_bound = max_theta_to_achieve_power( thetagrid( lower_theta_indices) ,...
                                          rejectiongrid( lower_theta_indices ),...
                                          power);

upper_bound = min_theta_to_achieve_power( thetagrid( upper_theta_indices) ,...
                                          rejectiongrid( upper_theta_indices ),...
                                          power);
 

end


function theta_min = min_theta_to_achieve_power( thetagrid, rejectiongrid, power)


has_power_indices = find( rejectiongrid >= power);


%If nothing meets the requisite power, return inf
if( isempty(has_power_indices) )
    theta_min = Inf;
    return;
end

%Take the minimum theta with at least the requisite power and the max theta
%without the requisite power (assuming that the min theta is not hte
%minimum in the grid -- in this case, return min of the grid). Do a linear
%interpolation between these two points
min_theta_above =  thetagrid( min(has_power_indices) );

if( min(has_power_indices) == 1)
    theta_min = min_theta_above;
else
    max_theta_below = thetagrid( min(has_power_indices) - 1);
    power_below = rejectiongrid( min(has_power_indices) -1 );
    power_above = rejectiongrid( min(has_power_indices) );
    
    theta_min = x_interpolated(max_theta_below, min_theta_above,power_below, power_above, power);


end

end

function theta_max = max_theta_to_achieve_power( thetagrid, rejectiongrid, power)


has_power_indices = find( rejectiongrid >= power);


%If nothing meets the requisite power, return -inf
if( isempty(has_power_indices) )
    theta_max = -Inf;
    return;
end

%Take the max theta with at least the requisite power and the min theta
%without the requisite power (assuming that the max theta is not the
%max in the grid -- in this case, return max of the grid). Do a linear
%interpolation between these two points
max_theta_above =  thetagrid( max(has_power_indices) );

if( max(has_power_indices) == 1)
    theta_max = max_theta_above;
else
    min_theta_below = thetagrid( max(has_power_indices) + 1);
    power_below = rejectiongrid( max(has_power_indices) + 1 );
    power_above = rejectiongrid( max(has_power_indices) );
    
    theta_max = x_interpolated(max_theta_above, min_theta_below,power_above, power_below, power);


end

end