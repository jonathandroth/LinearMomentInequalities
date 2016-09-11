

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
T = 27*1000 + burnout;

%g_vec = [31.818; 12.7; 21.495; 31.494; 43.462; 51.616; 54.277]/10;
g_vec = linspace(12.7, 54.277, 22)'/10; %Divide by 10 to put into 10k lb units

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

markup = .35;
avg_price = mean( [64510, 68012, 67644]);
observed_variance_quantity = var( [467; 502; 494]/1000);%divide by 1000 to get in millions


%Create a table of the target market shares

%latex2( '%.2f', {'Firm', 'Target Market Share', 'Target Avg. Products'}, f_names, f_share, avg_products_per_firm)
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

rng(0);
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






%% Calibrate sigma_nu and sigma_eps to match variance in 9-year averages
% 
% 
% sigma_nu_vec = [10;20;30;40];
% %sigma_eps_vec = sigma_nu_vec;
% 
% %sigma_nu_vec = 16;
% %sigma_eps_vec = 10;
% %sigma_nu_vec =  [10;11;12;13;14];
% %sigma_nu_vec = 10 + (0:10)'/10;
% sigma_eps_vec = [10;20;30;40];
% 
% %sigma_nu_vec = [0.01; 0.1; 0.5; 0.7; 1];
% %sigma_eps_vec = [.001];
% var_vec = NaN(size(sigma_nu_vec) );
% p_e = NaN(size(sigma_nu_vec,1), size(sigma_eps_vec,1));
% p_d = NaN(size(p_e));
% 
% for s = 1:length(sigma_nu_vec)
% 
%      
% parfor sprime = 1:length(sigma_eps_vec)
% 
%     %pct_done = (s-1)*sprime / (length(sigma_nu_vec) * length(sigma_eps_vec) )
%    
%  
% [mu_f, J_t, J_tm1] = mu_f_optimal( avg_products_per_firm, mu_f_vec, sigma_nu_vec(s), ...
%     sigma_eps_vec(sprime), zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, G_array);
% 
% var_vec(s,sprime) = nineyr_variance(J_t,burnout);
% 
% p_e(s,sprime) = mean( mean( mean( J_t(J_tm1 == 0),3 ), 2), 1 );
% p_d(s, sprime) = 1 - mean( mean( mean( J_t(J_tm1 == 1),3 ), 2), 1 );
% 
% end
% end
% var_vec
% %[p_e, p_d ]
% %p_e ./ (p_e + p_d)
% 
% %Choose the parameters that make the variance closest to
% %observed_3pd_variance
% 
% % [~,sigma_nu_ind] = min( min( abs(var_vec - observed_threepd_variance),[], 2) ) ;
% % [~,sigma_eps_ind] = min( min( abs(var_vec - observed_threepd_variance),[], 1) ) ;
% % 
% % sigma_nu = sigma_nu_vec(sigma_nu_ind);
% % sigma_eps = sigma_eps_vec(sigma_eps_ind);
% 


%% Calibrate sigma_zeta to match variance in pi

sigma_nu = 30;
sigma_eps = 30;

 [mu_f, J_t, J_tm1, Pi_star_array] = mu_f_optimal( avg_products_per_firm, mu_f_vec, sigma_nu, ...
     sigma_eps, zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, G_array);

%This function returns the pi_star array as a function of sigma_z. It sets
%sigma_zetaj= sigma_zetajft, and uses the shocks and the pi_star already
%calculated
pi_star_sigma_fn = @(sigma_z) pi_array_fn( Pi_star_array, Zetaj_shocks_array, Zetajft_shocks_array, sigma_z, sigma_z);

%This function takes the 9-year variance of the sum of the pi-stars for the
%products that were in the market
pi_star_variance_fn = @(sigma_z) nineyr_variance( J_t .* pi_star_sigma_fn(sigma_z) , burnout )  ;

target_variance_pi = markup^2 * avg_price^2 * observed_variance_quantity;

%Variance if there were no sigma_z as a fraction of target (check that htis is < 1)
nineyr_variance( Pi_star_array .* J_t, burnout) / target_variance_pi

myopts = optimset('TolX', 10^(-6), 'TolFun', 10^(-6) );
sigma_zeta = fzero( @(sigma_z) pi_star_variance_fn(sigma_z) - target_variance_pi, [0,target_variance_pi])


%% Make tables of calibrated params

%latex2( '%.2f', {'Firm', '$\mu_f$'}, f_names, mu_f)
%latex2( '%.2f', {'Parameter', 'Value'}, {'$\sigma_\nu$' ; '$\sigma_\epsilon$';'$\sigma_\zeta$' }, [sigma_nu; sigma_eps; sigma_zeta ])

%Save the calibrated parameters
save('../../Output/Simulated_Data/calibrated_params', 'sigma_nu', 'sigma_eps', 'sigma_zeta');

%% Simulate Data using Calibrated Parameters
% Now, we simulate one long chain

load('../../Output/Simulated_Data/calibrated_params');


%Re-set T
T = 50000 + burnout;

   
    %Reset F and G
    F_array = repmat( (1:F)',1, J, T);
    G_array = repmat( g_vec',F,1,T);

    % Re-draw the shocks for this simulation
    
        Eta_shocks_array = NaN(F,J,T);
        Epsilon_shocks_array = NaN(F,J,T);
        Zetaj_shocks_array = NaN(F,J,T);
        Zetajft_shocks_array = NaN(F,J,T);

        rng(1);
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
    
     
    % Draw product offerings and pi_star
        [J_t_array, J_tminus1_array, Pi_star_array] = calculate_offerings( sigma_nu, ...
        sigma_eps, zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, mu_f, G_array );
    
    % Construct pi from pi_star
        Pi_array = pi_array_fn( Pi_star_array, Zetaj_shocks_array, Zetajft_shocks_array, sigma_zeta, sigma_zeta);
    
    %Take burnout
        take_burnout_fn = @(X) X(:,:,(burnout+1):(size(X,3)) );
        
        J_t_array = take_burnout_fn(J_t_array);
        J_tminus1_array = take_burnout_fn(J_tminus1_array);
        Pi_star_array = take_burnout_fn(Pi_star_array);
        Pi_array = take_burnout_fn(Pi_array);
        Eta_shocks_array = take_burnout_fn( Eta_shocks_array);
        
        F_array = repmat( (1:F)',1, J, T-burnout);
        G_array = repmat( g_vec',F,1,T-burnout);
    
    
    % Save the desired values
        ds_name = '../../Output/Simulated_Data/Calibrated_SigmaZeta/ds_long';
        save( ds_name, 'J_t_array', 'J_tminus1_array', 'Pi_array', 'F_array', 'G_array', 'Pi_star_array', 'Eta_shocks_array');

    
    

%% Simulate data with various sigma_zeta/4

load('../../Output/Simulated_Data/calibrated_params');

sigma_zeta = sigma_zeta/4;

   
    %Reset F and G
    F_array = repmat( (1:F)',1, J, T);
    G_array = repmat( g_vec',F,1,T);

    % Re-draw the shocks for this simulation
    
        Eta_shocks_array = NaN(F,J,T);
        Epsilon_shocks_array = NaN(F,J,T);
        Zetaj_shocks_array = NaN(F,J,T);
        Zetajft_shocks_array = NaN(F,J,T);

        rng(1);
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
    
     
    % Draw product offerings and pi_star
        [J_t_array, J_tminus1_array, Pi_star_array] = calculate_offerings( sigma_nu, ...
        sigma_eps, zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, mu_f, G_array );
    
    % Construct pi from pi_star
        Pi_array = pi_array_fn( Pi_star_array, Zetaj_shocks_array, Zetajft_shocks_array, sigma_zeta, sigma_zeta);
    
    %Take burnout
        take_burnout_fn = @(X) X(:,:,(burnout+1):(size(X,3)) );
        
        J_t_array = take_burnout_fn(J_t_array);
        J_tminus1_array = take_burnout_fn(J_tminus1_array);
        Pi_star_array = take_burnout_fn(Pi_star_array);
        Pi_array = take_burnout_fn(Pi_array);
        Eta_shocks_array = take_burnout_fn( Eta_shocks_array);
        
        F_array = repmat( (1:F)',1, J, T-burnout);
        G_array = repmat( g_vec',F,1,T-burnout);
    
    
    % Save the desired values
        ds_name = strcat( '../../Output/Simulated_Data/Calibrated_SigmaZeta_Over4/ds_long');
        save( ds_name, 'J_t_array', 'J_tminus1_array', 'Pi_array', 'F_array', 'G_array', 'Pi_star_array', 'Eta_shocks_array');

    
   
display('Simulation Complete.');

%% Simulate data with sigma_zeta = 0


load('../../Output/Simulated_Data/calibrated_params');

sigma_zeta = 0;

   
    %Reset F and G
    F_array = repmat( (1:F)',1, J, T);
    G_array = repmat( g_vec',F,1,T);

    % Re-draw the shocks for this simulation
    
        Eta_shocks_array = NaN(F,J,T);
        Epsilon_shocks_array = NaN(F,J,T);
        Zetaj_shocks_array = NaN(F,J,T);
        Zetajft_shocks_array = NaN(F,J,T);

        rng(1);
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
    
     
    % Draw product offerings and pi_star
        [J_t_array, J_tminus1_array, Pi_star_array] = calculate_offerings( sigma_nu, ...
        sigma_eps, zeros(F,J) , Epsilon_shocks_array, Eta_shocks_array, mu_f, G_array );
    
    % Construct pi from pi_star
        Pi_array = pi_array_fn( Pi_star_array, Zetaj_shocks_array, Zetajft_shocks_array, sigma_zeta, sigma_zeta);
    
    %Take burnout
        take_burnout_fn = @(X) X(:,:,(burnout+1):(size(X,3)) );
        
        J_t_array = take_burnout_fn(J_t_array);
        J_tminus1_array = take_burnout_fn(J_tminus1_array);
        Pi_star_array = take_burnout_fn(Pi_star_array);
        Pi_array = take_burnout_fn(Pi_array);
        Eta_shocks_array = take_burnout_fn( Eta_shocks_array);
        
        F_array = repmat( (1:F)',1, J, T-burnout);
        G_array = repmat( g_vec',F,1,T-burnout);
    
    
    % Save the desired values
        ds_name = strcat( '../../Output/Simulated_Data/SigmaZeta_Equal0/ds_long');
        save( ds_name, 'J_t_array', 'J_tminus1_array', 'Pi_array', 'F_array', 'G_array', 'Pi_star_array', 'Eta_shocks_array');

    
    

display('Simulation Complete.');