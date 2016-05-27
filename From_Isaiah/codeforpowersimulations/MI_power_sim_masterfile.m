%parpool('local',24)

mu_other_grid=[0 -1 -2 -3 -5 -7];
rho_grid=0;
k_grid=[2 5 10 20 100]; %5 states

%Set grid of alternatives considered
%mu_grid=-10:.1:10;
mu_grid=-10:.5:10;

%set number of simulations
sims=10^3;
crit_sims=10^4;

power_store_huge=zeros(length(mu_grid),length(mu_grid),length(mu_other_grid),...
    length(rho_grid),length(k_grid));
Twostep_power_store_huge=zeros(length(mu_grid),length(mu_grid),length(mu_other_grid),...
    length(rho_grid),length(k_grid));
LF_power_store_huge=zeros(length(mu_grid),length(mu_grid),length(mu_other_grid),...
    length(rho_grid),length(k_grid));
RSW_power_store_huge=zeros(length(mu_grid),length(mu_grid),length(mu_other_grid),...
    length(rho_grid),length(k_grid));

for l3=1:length(k_grid)
    for l2=1:length(rho_grid)
        if l3==1;
            l1=1;
            mu_other=mu_other_grid(1);
            l1
            l2
            l3
            k=k_grid(l3);
            rho=rho_grid(l2);
            
            MI_power_sim_subfile
            
            power_store_huge(:,:,l1,l2,l3)=power_store;
            Twostep_power_store_huge(:,:,l1,l2,l3)=Twostep_power_store;
            LF_power_store_huge(:,:,l1,l2,l3)=LF_power_store;
            RSW_power_store_huge(:,:,l1,l2,l3)=RSW_power_store;       
            save MI_power_sim_temp_backup1
        else
            for l1=1:length(mu_other_grid)
                
                l1
                l2
                l3
                mu_other=mu_other_grid(l1);
                k=k_grid(l3);
                rho=rho_grid(l2);
                
                MI_power_sim_subfile
                
                power_store_huge(:,:,l1,l2,l3)=power_store;
                Twostep_power_store_huge(:,:,l1,l2,l3)=Twostep_power_store;
                LF_power_store_huge(:,:,l1,l2,l3)=LF_power_store;
                RSW_power_store_huge(:,:,l1,l2,l3)=RSW_power_store;
                
                save MI_power_sim_temp_backup1
            end
        end        
        save MI_power_sim_temp_backup2
    end
    save MI_power_sim_temp_backup3
end
