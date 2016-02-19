% This is the primary function for simulating a draw of the data (i.e. all
% T firms, G weights, and the corresponding delta pi's)

%It first creates the J distribution for each firm, by calling the
%sim_offering_firmf function. Then it draws the delta pi's.

function [t,f, g, C_mat, Delta_pi_mat] = sim_data()

global theta_c lambda theta_g T G F p_d p_o_f gamma; 





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

end 
