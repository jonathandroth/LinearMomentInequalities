function [mu_f,J,J_tminus1, Pi_star] = mu_f_optimal( avg_products_per_firm, mu_f_guess, sigma_nu, ...
    sigma_epsilon, J_0 , Epsilon_shocks_array, Eta_shocks_array, Eta_t_array,  G_array)

global burnout;
global F;
global theta_c;

%Create a function that computes the J_array for a given value of mu_f,
%holding constant all the other parameters
J_array_muf = @(mu_f) calculate_offerings( sigma_nu, ...
    sigma_epsilon, J_0 , Epsilon_shocks_array, Eta_shocks_array, Eta_t_array,  mu_f, G_array);

%Create a function that computes the mean products by firm using the output
%from the J_array_muf function
mean_by_firm_muf = @(mu_f) mean_products_byfirm( J_array_muf(mu_f) , burnout );


%firstelement = @(x) x(1);
%mu1_obj = @(mu_f1) firstelement( mean_by_firm_muf( [mu_f1; zeros(8,1) ]) )- avg_products_per_firm(1) ;

%Solve for the value of mu_f that sets the mean_by_firm in the data equal
%to the desired avg_products_per_firm

%myopts = optimset('TolX', 10^(-4), 'TolFun', 10^(-4) );
%mu_f = fzero( mu1_obj, [0;200], myopts);
%mu_f = binary_search( mu1_obj, [0;200]);

%myopts = optimset('TolX', 10^(-4), 'TolFun', 10^(-4), 'FinDiffRelStep', .5 );
%mu_f = fsolve( @(mu_f) mean_by_firm_muf(mu_f) - avg_products_per_firm, repmat(100,9,1), myopts);

mu_f = NaN(F,1);

theta_c_local = theta_c;
for f = 1:F
   
    obj_f = @(mu_f) mean_products_firm_f( f, mu_f, mean_by_firm_muf ) - avg_products_per_firm(f);
    myopts = optimset('TolX', 10^(-4), 'TolFun', 10^(-4) );
    mu_f(f) = fzero( obj_f, [theta_c_local - 2*sigma_nu - 2*sigma_epsilon ;theta_c_local + 2*sigma_nu + 2*sigma_epsilon],myopts );
    
end

[J, J_tminus1, Pi_star] = J_array_muf(mu_f);

end


%This helper function returns the average products for firm f given a mean
%mu_f for firm f. It takes as arguments the firm number f, the mean mu_f,
%and the function mean_by_firm_muf, which takes a mean for all of the firms
%and returns the market share for all of the firms
function firmf_mean = mean_products_firm_f( f, mu_f, mean_by_firm_muf )
    global F;

    e_f = zeros(F,1);
    e_f(f) = mu_f;
    
    all_firm_means = mean_by_firm_muf(e_f);
    
    firmf_mean = all_firm_means(f);
    
end