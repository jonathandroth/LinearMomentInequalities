%%This function simulates the offerings for firm firmnum

%It takes in the transition probabilities p_o and p_d, and the firmnumber
%firmnum

%It returns data in long format. The vectors t, f and g indicate the time, firmnumber, and
%weight class of the product. The vectors in_t and in_tminus1 indicate
%whether the truck with weight g was offered in period t and period t-1

function [t,f, g, in_t , in_tminus1] = sim_offerings_firmf(firmnum, p_o, p_d)

global T;
global G;


%Create a matrix of size (T+1) x G
%%Each row of the matrix will represent a time period
%%Each column of the matrix will represent a weight class.
%%Each entry will indicate whether the firm offered that weight in the given time period
Jft = NaN(T+1, G);

%Draw the first row of the matrix (for t=0) from the stationary
%distribution given p_o and p_d
Jf0 =  rand(1, G) <( (p_o / (p_o + p_d)) );
Jft(1, :) = Jf0;

%For period 1,...,t; simulate the next period's offerings using the
%previous period's offerings and the transition probs (uses the helper fn
%next_Jft)

for t= 1:T

    Jft(t+1,:) = next_Jft( Jft(t,:) , p_o, p_d );

end

%%%%% Reshape data %%%%%%%
% The following block of code reshapes the data into a "long format." That
% is, I will have one column indicating the period t, one column indicating
% hte weight class g. There will then be columns in_t and in_tminus1,
% indicating whether the truck with weight g was in in period t and whether
% it was in in period t-1


%Create two matrices, one with all the observations from period t=1 through
%T, and one with all the observations from period t=0 through t-1

%These will be reshaped into vectors in_t and in_tminus1

Jft_thisperiod = Jft( 2:(T+1), :);
Jft_lastperiod = Jft( 1:T, :);

%Create a matrix indicating the year of each observation in Jft_thisperiod
t_mat = repmat( (1:T)', 1, G); 
%Create a matrix indicating the g of each observation in Jft_thisperiod
g_mat = repmat( (1:G), T, 1);

%Reshape all of the above matrices into vectors
in_t = Jft_thisperiod(:);
in_tminus1 = Jft_lastperiod(:);

t = t_mat(:);
g = g_mat(:);

f = repmat(firmnum, size(t) );

end
