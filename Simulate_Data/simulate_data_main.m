

%% Set the parameter values
global theta_c;
theta_c = 129.73;

global lambda;
lambda = 0.386; %should this be negative? Should this be one-over this? %Think we are good.
%lambda = 1;

global theta_g;
theta_g = -21.38;


global burnout;
burnout = 1000;


global T;
T = 27*100 + burnout;

%g_vec = [31.818; 12.7; 21.495; 31.494; 43.462; 51.616; 54.277];
g_vec = linspace(12.7, 54.277, 12)';


global J;
J = size(g_vec,1);
 
global F;
F = 9;


f_names = {'Chrysler'; 'Ford' ;'Daimler'; 'GM'; 'Hino'; 'International'; 'Isuzu'; ...
    'Paccar'; 'Volvo'} ;

f_share = [5.56; 12.5; 12.5+4.17 + 1.39; 6.94; 4.17; 16.67; 6.94; 6.94+11.11; 2.78...
    + 4.17+4.17]/100; 




mu_f_vec = [repmat(100,F-1,1); 50];


global avg_total_yearly_products;
avg_total_yearly_products = 48;

avg_products_per_firm = f_share * avg_total_yearly_products;

observed_9yr_totals = [8 + 21.3 + 8.3; 7.9 + 24.8 + 10.8; 8.3 + 27.1 + 9.6];
observed_threepd_variance = var( observed_9yr_totals);
%% Make covariate arrays

%Make arrays corresponding to the covariates for the products in the shocks
%arrays. Again, the arrays are F x J x T

F_array = repmat( (1:F)',1, J, T);
G_array = repmat( g_vec',F,1,T);





%% Draw the shocks


%The shocks will be stored in a three dimensional array of size F X J x T
% The first dimension (rows) will represent different firms, the second
% dimension (columns) will represent the different products, and the third
% dimension (depth) will represent different time periods

%I will take the draws of the shocks from standard normals. To go from
% Eta_shocks to the Eta_draws, I will multiply by sigma_eta. Likewise, to
% go from the Epsilon_shocks to the Epsilon_draws,  I will multiply by sigma_f
% and then add mu_f to each row. I start with the standard shocks so that
% we can use the same draws regardless of the parameters


Eta_shocks_array = NaN(F,J,T);
Epsilon_shocks_array = NaN(F,J,T);
Zetaj_shocks_array = NaN(F,J,T);
Zetajft_shocks_array = NaN(F,J,T);

%rng(1234);
for t= 1:T
    %For eta, take one draw for each j for period t
    %Right now, assuming that the shocks are uncorrelated across J. Could
    %relax this (change away from eye(J) )
    eta_shocks_t = mvnrnd( zeros(1,J), eye(J) ); 
    
    %Do the same for zeta_j
    zetaj_shocks_t = mvnrnd( zeros(1,J), eye(J) );
    
    %Create and store the eta draw matrix for period t.
    %The eta draw is the same for product j for each firm, so the eta_shocks
    %matrix for period t has constant columns (this is done in the repmat)
    Eta_shocks_array(:,:,t) = repmat( eta_shocks_t, F, 1);
   
    %Do the same for zeta_j
    Zetaj_shocks_array(:,:,t) = repmat( zetaj_shocks_t, F, 1);
    
    %Draw the epsilon_tfj. These are assumed to be independent across all
    %draws
    Epsilon_shocks_array(:,:,t) = randn(F,J);
    
    %Do the same for zeta_jft
    Zetajft_shocks_array(:,:,t) = randn(F,J);
   
end




%% Draw data
sigma_nu = 50;
sigma_zetaj = 50;
% 
% tic;
% [J_t_startat1, J_tm1_startat1] = calculate_offerings( sigma_nu, ...
%     sigma_epsilon, ones(F,J) , Epsilon_shocks_array, Eta_shocks_array, Zetaj_shocks_array,Zetajft_shocks_array, mu_f_vec, G_array, F_array );
% toc;
% 
% 
% [J_t_startat0, J_tm1_startat0] = calculate_offerings( sigma_nu, ...
%     sigma_epsilon, zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, Zetaj_shocks_array,Zetajft_shocks_array, mu_f_vec, G_array, F_array );
% 

sigma_nu_vec = [100;250;500;750;1000];
sigma_eps_vec = [100;200;300];
%sigma_eps_vec = [100;500;1000];
var_vec = NaN(size(sigma_nu_vec) );
p_e = NaN(size(sigma_nu_vec,1), size(sigma_eps_vec,1));
p_d = NaN(size(p_e));

for s = 1:length(sigma_nu_vec)
for sprime = 1:length(sigma_eps_vec)

[mu_f, J_t, J_tm1] = mu_f_optimal( avg_products_per_firm, mu_f_vec, sigma_nu_vec(s), ...
    sigma_eps_vec(sprime), zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, G_array);

var_vec(s,sprime) = nineyr_variance(J_t,burnout);

p_e(s,sprime) = mean( mean( mean( J_t(J_tm1 == 0),3 ), 2), 1 );
p_d(s, sprime) = 1 - mean( mean( mean( J_t(J_tm1 == 1),3 ), 2), 1 );

end
end
var_vec
%[p_e, p_d ]
%p_e ./ (p_e + p_d)

%Choose the parameters that make the variance closest to
%observed_3pd_variance

[~,sigma_nu_ind] = min( min( abs(var_vec - observed_threepd_variance),[], 2) ) ;
[~,sigma_eps_ind] = min( min( abs(var_vec - observed_threepd_variance),[], 1) ) ;

sigma_nu = sigma_nu_vec(sigma_nu_ind);
sigma_eps_ = sigma_eps_vec(sigma_eps_ind);



%%
