%This script should be called to combine the moments related to Theta_g for
%A and X_T and Y_T. This is to be called both when creating hte moments in
%the main script (lp_create_moments_and_covariances) and when finding the
%identified set in lp_compute_confidence_sets

        moment_nums = 1:size(A,2);
        moment_nums_mod6 = mod2( moment_nums,6 );
        theta_g_cols = moment_nums_mod6 == 5 | moment_nums_mod6 == 6;       
        moment_nums_theta_g_cols = mod2( moment_nums( theta_g_cols), size(A,2) / num_F_groups);
        moment_nums( theta_g_cols) = moment_nums_theta_g_cols;
         
        A = grpstats2( A' , moment_nums')';
        
        moment_nums = 1:size(X_T,1);
        moment_nums_mod6 = mod2( moment_nums,6 );
        theta_g_cols = moment_nums_mod6 == 5 | moment_nums_mod6 == 6;       
        moment_nums_theta_g_cols = mod2( moment_nums( theta_g_cols), size(X_T,1) / num_F_groups);
        moment_nums( theta_g_cols) = moment_nums_theta_g_cols;
        
        X_T = grpstats2( X_T , moment_nums');
        Y_T = grpstats2( Y_T , moment_nums');
        Y_wide = grpstats2( Y_wide', moment_nums')'; 
