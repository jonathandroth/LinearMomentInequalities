

%% Set the parameter values
global theta_0;
theta_0 = 129.73;

global lambda;
lambda = 0.386; %should this be negative? Should this be one-over this?

global theta_g;
theta_g = 21.38;



global T;
T = 10000;

global G;
G = 10;

%% Simulate the J distribution

Jft = NaN(T+1, G);

p_o = .5;
p_d = .25;

Jf0 =  rand(1, G) <( (p_o / (p_o + p_d)) );

Jft(1, :) = Jf0;
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


%% Test this
mean( in_t) 
mean( in_t( in_tminus1 == 1) )
mean( in_t( in_tminus1 == 0) )

%% Draw Delta_Pi | J





%% Construct the sample moments