function [CS_test,QLR_stat, CS_crit] = CoxAndShi_NoNuisance(Y, A, b, Sigma, alpha)
        
        %This function takes a vector Y, assumed to be ~normally
        %distributed with variance Sigma, and runs the CS test for the null
        %hypothesis that H_0: A E[Y] <= b

        Sigmahat_inv = Sigma^-1;
        options=optimset('Display','off','TolX',1e-7);
            %Compute the QLR_stat and optimizer
            % Note that here we use a re-parameterization of the
            % optimization, where the optimization variable is X = (Y-mu)
            
        [t_hat, QLR_stat]=quadprog((Sigmahat_inv+Sigmahat_inv'),[],-A,b-A*Y,[],[],[],[],[],options);
        %t_hat_renorm=Y-t_hat;
        moments = A*(Y-t_hat) - b;
        %binding_count=max(sum(moments>=-10^-6),1);
        binding_index=find(moments>=-10^-6);
        binding_count = max( rank(A(binding_index,:)) , 1);
        CS_crit=chi2inv(1-alpha,binding_count);

        %Comment out refinement and just do the unrefined test for now
%         if binding_count~=1
%             CSR_crit=CS_crit;
%         else
%             binding_ID=(t_hat_renorm>=-10^-4);
%             binding_norm=sqrt(binding_ID'*Sigmahat*binding_ID);
%             num=binding_norm*(-Y);
%             denom=binding_norm*diag(Sigmahat).^0.5-Sigmahat*binding_ID;
%             rat=num./denom;
%             rat(binding_ID==1)=[];
%             tauhat=min(rat);
%             betahat=2*alpha*normcdf(tauhat);
%             CSR_crit=chi2inv(1-betahat,binding_count);
%         end
        
        CS_test=QLR_stat>CS_crit;
        %CSR_test=QLR_stat>CSR_crit; 
        
end        