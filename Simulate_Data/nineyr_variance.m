function [mean_threeperiod_variance] = nineyr_variance( J_t, burnout )
    
%Take the array of product offerings excluding the burnout
J = J_t(:,:, (burnout+1):(size(J_t,3)) );

%Collapse to a vector of total products offered per year
N_t = sum( sum(J,1) ,2 );
N_t = N_t(:);
%Create an index vector for the 9 year periods (i.e. a vector that is 1
%9-times, 2 9-times, etc.

nineyear_index = repmat( (1:size(N_t,1)/9 ), 9,1);
nineyear_index = nineyear_index(:);

%Take the mean in the 9 year periods
nineyear_means = accumarray( nineyear_index, N_t, [], @mean);

%Create an index indicating which group of 3 9-year periods an observation
%belongs to (i.e. whether it represents the mean from the first 27 years,
%the second, etc). Practically, this vector is 1 3-times, 2 3-times, etc.
threeperiod_index = repmat( 1:( size(nineyear_means,1)/3), 3,1);
threeperiod_index = threeperiod_index(:);

%Compute the variance of the 9-year totals for each 27-year period
threeperiod_variances = accumarray( threeperiod_index, nineyear_means, [], @var );

%Compute the mean of these variances
mean_threeperiod_variance = mean(threeperiod_variances);
end