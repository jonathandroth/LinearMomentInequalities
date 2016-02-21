

%% Set the parameter values
global theta_c;
theta_c = 129.73;

global lambda;
lambda = 0.386; %should this be negative? Should this be one-over this?

global theta_g;
theta_g = 21.38;


global burnout;
burnout = 1000;


global T;
T = 27 + burnout;

global J;
J = 10;
 
global F;
F = 9;


f_names = {'Chrysler'; 'Ford' ;'Daimler'; 'GM'; 'Hino'; 'International'; 'Isuzu'; ...
    'Paccar'; 'Volvo'} ;

f_share = [5.56; 12.5; 12.5+4.17 + 1.39; 6.94; 4.17; 16.67; 6.94; 6.94+11.11; 2.78...
    + 4.17+4.17]/100; 



mu_f_vec = repmat(100,F,1);

%% Make covariate arrays

%Make arrays corresponding to the covariates for the products in the shocks
%arrays. Again, the arrays are F x J x T

F_array = repmat( (1:F)',1, J, T);
G_array = repmat( (1:J),F,1,T);





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


Mu_f_array = repmat(mu_f_vec, 1,J,T);


%% Draw data
sigma_nu = 50;
sigma_epsilon = 50;


[J_t_startat1, J_tm1_startat1] = calculate_offerings( sigma_nu, ...
    sigma_epsilon, ones(F,J) , Epsilon_shocks_array, Eta_shocks_array, Zetaj_shocks_array,Zetajft_shocks_array, Mu_f_array, G_array, F_array );

[J_t_startat0, J_tm1_startat0] = calculate_offerings( sigma_nu, ...
    sigma_epsilon, zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, Zetaj_shocks_array,Zetajft_shocks_array, Mu_f_array, G_array, F_array );


mean(J_t_startat1(:))
mean(J_t_startat0(:))


meanyrt = @(J,t) mean( mean(J(:,:,t)) );

meanyr_J0 = @(t) meanyrt(J_t_startat0, t);
meanyr_J1 = @(t) meanyrt(J_t_startat1, t);

meanyr_J0_tm1 = @(t) meanyrt(J_tm1_startat0, t);
meanyr_J1_tm1 = @(t) meanyrt(J_tm1_startat1, t);


J0_means = arrayfun( meanyr_J0 , (1:T)' );
J1_means = arrayfun( meanyr_J1 , (1:T)' );

J0_means_tm1 = arrayfun( meanyr_J0_tm1 , (1:T)' );
J1_means_tm1 = arrayfun( meanyr_J1_tm1 , (1:T)' );


plot( [(1:T)', (1:T)'], [J0_means, J1_means] )

[J0_means_tm1, J1_means_tm1]
%plot( [(1:T)'], [J0_means] )


mean( mean(J_t_startat0, 3), 1  ) 