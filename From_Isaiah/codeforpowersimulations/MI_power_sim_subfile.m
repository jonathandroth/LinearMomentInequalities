Sigma=ones(k,k)*rho+eye(k)*(1-rho);
Sigma_rt=Sigma^.5;
Sigma_diag=diag(Sigma);

alpha=0.05;
beta=0.005;
two_step_alpha=alpha-beta;

RSW_beta=.005;
RSW_alpha=alpha-RSW_beta;

noise_draw_crit=Sigma_rt*randn(k,crit_sims);
LF_crit=quantile(max(noise_draw_crit,[],1),1-alpha);
Pretest_crit=quantile(max(noise_draw_crit,[],1),1-beta);
RSW_pretest_crit=quantile(max(noise_draw_crit,[],1),1-RSW_beta);

noise_draw=Sigma_rt*randn(k,sims);
power_store=zeros(length(mu_grid),length(mu_grid));
LF_power_store=zeros(length(mu_grid),length(mu_grid));
Twostep_power_store=zeros(length(mu_grid),length(mu_grid));
RSW_power_store=zeros(length(mu_grid),length(mu_grid));
for n1=1:length(mu_grid)
    for n2=1:length(mu_grid)
       mu=zeros(k,1);
       mu(1)=mu_grid(n1);
       mu(2)=mu_grid(n2);
       if k>2
          mu(3:end)=mu_other;
       end
       
       draw=noise_draw+repmat(mu,1,sims);
       
       test=zeros(sims,1);
       Twostep_test=zeros(sims,1);
       LF_test=zeros(sims,1);
       RSW_test=zeros(sims,1);
       for m=1:sims
       %parfor m=1:sims
         
       %Draw vector of moments
       X=draw(:,m);
       %calcualte max of moments
       max_stat=max(X);
       %pick out moment where max occurs
       Max_ind=X==max_stat;
       %Calculate conditioning statistics
       xi=(X(Max_ind==0)-Sigma(Max_ind==0,Max_ind==1)./Sigma_diag(Max_ind==0)*max_stat)...
           ./(1-Sigma(Max_ind==0,Max_ind==1)./Sigma_diag(Max_ind==0));
       J=Sigma(Max_ind==0,Max_ind==1)./Sigma_diag(Max_ind==0);
       V_lo=max(xi(J<1));
       V_up=min(xi(J>1));
       if max(J)<=1
          V_up=inf; 
       end
       %Calculate condtional p value
       T=1-(normcdf(max_stat)-normcdf(V_lo))./...
           (normcdf(V_up)-normcdf(V_lo));
       
       %Since running into numerical precision issues in calculating
       %conditinal critical value, also calculate upper bound on p value
       T_bound=1-(((4+V_lo.^2).^.5+V_lo).^-1-((2+max_stat.^2).^.5+max_stat).^-1.*...
           exp(-.5*(max_stat.^2-V_lo.^2)))./(((2+V_lo.^2).^.5+V_lo).^-1-...
           ((4+V_up.^2).^.5+V_up).^-1.*exp(-.5*(V_up.^2-V_lo.^2)));
       T_bound(min(V_lo,max_stat)<0)=1;
       
       T=min(T,T_bound);
       
       test(m,1)=T<alpha;
       Twostep_test(m,1)=max(T<two_step_alpha,max_stat>Pretest_crit);
       LF_test(m,1)=max_stat>LF_crit;
       
       %Calculate Romano Shaikh and Wolf Test
       lambda_RSW=X+RSW_pretest_crit;
       lambda_RSW=min(lambda_RSW,zeros(k,1));
       RSW_noise=noise_draw_crit+repmat(lambda_RSW,1,crit_sims);
       RSW_crit=quantile(max(RSW_noise,[],1),1-RSW_alpha);
       
       RSW_test(m,1)=(max_stat>RSW_crit)*max(lambda_RSW==0);  
       end
       power_store(n1,n2)=mean(test);   
       Twostep_power_store(n1,n2)=mean(Twostep_test);
       LF_power_store(n1,n2)=mean(LF_test);
       RSW_power_store(n1,n2)=mean(RSW_test);
    end
end

power_diff=power_store-LF_power_store;