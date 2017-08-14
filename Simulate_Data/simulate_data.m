function [J_t_array, J_tminus1_array, Pi_array, F_array, G_array, Pi_star_array, Eta_jt_shocks_array, Eta_t_vec] = simulate_data( seed, F , J ,T, burnout, sigma_nu, sigma_eps, sigma_w, sigma_zetaj, sigma_zetajft, rho, lambda,theta_c,theta_g, g_vec, mu_f)
      
    G_array = repmat( g_vec',F,1,T);
     
     %G_array = repmat( 10*g_vec,F,1,T); %Make G larger than in the calibration

    % Re-draw the shocks for this simulation
    
        Eta_jt_shocks_array = NaN(F,J,T);
        Eta_t_vec = NaN(T,1);
        Epsilon_shocks_array = NaN(F,J,T);
        Zetaj_shocks_array = NaN(F,J,T);
        Zetajft_shocks_array = NaN(F,J,T);

        rng(seed);
        
        
        eta_tminus1 = 0;
        
        for t= 1:T
            %For eta, take one draw for each j for period t
            %Right now, assuming that the shocks are uncorrelated across J. Could
            %relax this (change away from eye(J) )
            eta_jt_shocks_t = mvnrnd( zeros(1,J), eye(J) ); 

            
            %Draw the innovation in eta_t, i.e. the market-wide time shock
            w_t = sigma_w * randn(1);
            
            %Update eta_t
            eta_t = rho * eta_tminus1 + w_t;
            eta_tminus1 = eta_t;
            
            %Draw zeta_j
            zetaj_shocks_t = mvnrnd( zeros(1,J), eye(J) );

            %Create and store the eta draw matrix for period t.
            %The eta draw is the same for product j for each firm, so the eta_shocks
            %matrix for period t has constant columns (this is done in the repmat)
            Eta_jt_shocks_array(:,:,t) = repmat( eta_jt_shocks_t, F, 1);

            %Do the same for zeta_j
            Zetaj_shocks_array(:,:,t) = repmat( zetaj_shocks_t, F, 1);

            %Draw the epsilon_tfj. These are assumed to be independent across all
            %draws
            Epsilon_shocks_array(:,:,t) = randn(F,J);

            %Do the same for zeta_jft
            Zetajft_shocks_array(:,:,t) = randn(F,J);
            
            %Store the eta_t shocks in a vector
            Eta_t_vec(t,1) = eta_t;
 
        end
    
     
    % Draw product offerings and pi_star
        [J_t_array, J_tminus1_array, Pi_star_array] = calculate_offerings( sigma_nu, ...
        sigma_eps, zeros(F,J) , Epsilon_shocks_array, Eta_jt_shocks_array, Eta_t_vec, mu_f, G_array,...
        T,F,J,lambda,theta_c,theta_g);
    
    
    % Construct pi from pi_star
        Pi_array = pi_array_fn( Pi_star_array, Zetaj_shocks_array, Zetajft_shocks_array, sigma_zetaj, sigma_zetajft);
    
    %Take burnout
        take_burnout_fn = @(X) X(:,:,(burnout+1):(size(X,3)) );
        
        J_t_array = take_burnout_fn(J_t_array);
        J_tminus1_array = take_burnout_fn(J_tminus1_array);
        Pi_star_array = take_burnout_fn(Pi_star_array);
        Pi_array = take_burnout_fn(Pi_array);
        Eta_jt_shocks_array = take_burnout_fn( Eta_jt_shocks_array);
        Eta_t_vec = Eta_t_vec( (burnout+1):end) ;
        
        F_array = repmat( (1:F),1, J, T-burnout);
        G_array = repmat( g_vec,F,1,T-burnout);

end