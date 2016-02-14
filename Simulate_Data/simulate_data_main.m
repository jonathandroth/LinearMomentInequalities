

%% Set the parameter values
global theta_c;
theta_c = 129.73;

global lambda;
lambda = 0.386; %should this be negative? Should this be one-over this?

global theta_g;
theta_g = 21.38;



global T;
T = 30;

global G;
G = 10;

global F;
F = 9;

f_names = {'Chrysler'; 'Ford' ;'Daimler'; 'GM'; 'Hino'; 'International'; 'Isuzu'; ...
    'Paccar'; 'Volvo'} ;

f_share = [5.56; 12.5; 12.5+4.17 + 1.39; 6.94; 4.17; 16.67; 6.94; 6.94+11.11; 2.78...
    + 4.17+4.17]/100; 

global p_d;
p_d = .25;

total_avg_products = 37;

% The number of products offered on average by each firm is G * p_o / (p_o
% + p_d). We this required that G * p_o / (p_o
% + p_d) = share * total_avg_products

% This implies p_o = share * total_avg_products * p_d / (G - share *
% total_avg_products) 

global p_o_f;
p_o_f = (f_share * total_avg_products * p_d) ./ (G - f_share * total_avg_products); 


global gamma;
gamma = [1;1;1;0;0;0];

%% Simulate the J distribution


%%Loop over all the firms and combine the simulated data for each firm (by appending the vectors);



[t,f, g, in_t , in_tminus1] = sim_offerings_firmf(1, p_o_f(1), p_d);

for firm = 2:F
    
    
    [t_f,f_f, g_f, in_t_f , in_tminus1_f] = sim_offerings_firmf(firm, p_o_f(firm), p_d); 
    
    t = [t;t_f];
    f = [f;f_f];
    g = [g;g_f];
    in_t = [in_t; in_t_f];
    in_tminus1 = [in_tminus1; in_tminus1_f];
    
end
    
        
%% Create C_ift

C_mat = NaN( size(t,1) , 6);

C_mat(:,1) = (in_t == 1 & in_tminus1 == 1);
C_mat(:,2) = (in_t == 0 & in_tminus1 == 1);
C_mat(:,3) = (in_t == 1 & in_tminus1 == 0);
C_mat(:,4) = (in_t == 0 & in_tminus1 == 0);
C_mat(:,5) = (in_t == 1 & in_tminus1 == 0);
C_mat(:,6) = (in_t == 1 & in_tminus1 == 0);



%Check the fraction of years that have at least one of each type of
%conditioning set
%mean( grpstats( C_mat, t, {'mean'}) > 0 )
% [mean, grps] = grpstats( C_mat, {t, g} , {'mean'})


%% Draw Delta_Pi | J
errors = mvnrnd(gamma , eye(6) , size(C_mat,1) );

Delta_pi_mat = NaN( size(errors) );

Delta_pi_mat(:,1) = errors(:,1) + lambda * theta_c + lambda * theta_g * g;
Delta_pi_mat(:,2) = errors(:,2) - lambda * theta_c -  lambda * theta_g * g;
Delta_pi_mat(:,3) = errors(:,3) + theta_c +  theta_g * g;
Delta_pi_mat(:,4) = errors(:,4) - theta_c -  theta_g * g;
Delta_pi_mat(:,5) = errors(:,5) + theta_g;
Delta_pi_mat(:,6) = errors(:,6) - theta_g;

%set delta pi's to NaN where the correspond C_i = 0
Delta_pi_mat( C_mat == 0) = NaN;










%% Construct the sample moments


moments = sample_moments(Delta_pi_mat,t, g, lambda, theta_c, theta_g);

