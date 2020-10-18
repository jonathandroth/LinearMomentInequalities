%This function calculates the average area of the accepted region given a
%2-d grid of rejection probabilities

%Its inputs are:

%rejection_grid: a 2d grid of rejection probabilities
% dim1, dim2: the dimensions of the space between the gridpoints in both
% directions

%Its output is

%accepted_area: the estimated average area of the accepted region. This is computed
%by summing the average acceptance probability at each point times the area
% of the box surrounding that point (i.e. dim1 * dim2)

function [accepted_area] = calculate_accepted_area( rejection_grid, dim1, dim2)



acceptance_prob = 1 - rejection_grid;

area = dim1 * dim2;

accepted_area = sum(acceptance_prob(:) ) * area;

end