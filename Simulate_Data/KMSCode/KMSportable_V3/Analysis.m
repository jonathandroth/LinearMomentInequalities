

%% Count feasible points
count_flagfeas=0;
ind_flagfeas = ones(Nmc,1);
for mm = 1:Nmc
    if KMS_output{mm}.flag_feas == 0
        count_flagfeas = count_flagfeas+1;
        ind_flagfeas(mm,1) = 0;
    end
end

count_flagLB=0;
count_flagUB = 0;
for mm = 1:Nmc
    if KMS_output{mm}.flag_feas == 1 
        if KMS_output{mm}.flagL_EAM==0
            count_flagLB = count_flagLB+1;
        end
        if KMS_output{mm}.flagU_EAM==0
            count_flagUB = count_flagUB+1;
        end
    end
end

%% KMS_CI (with feasible point)
med_KMS = median(KMS_CI(ind_flagfeas==1,:));

KMS_CI_feas = NaN(Nmc-count_flagfeas,2);
c_theta = NaN(Nmc-count_flagfeas,2);
totaltime = NaN(Nmc-count_flagfeas,1);
time_EAM    = NaN(Nmc-count_flagfeas,2);  
count = 1;
for mm = 1:Nmc
    if KMS_output{mm}.flag_feas == 1
        KMS_CI_feas(count,:) = KMS_CI(mm,:);
        c_theta(count,1) = KMS_output{mm}.cL_EAM;
        c_theta(count,2) = KMS_output{mm}.cU_EAM;
        time_EAM(count,1)= KMS_output{mm}.timeL_EAM;
        time_EAM(count,2)= KMS_output{mm}.timeU_EAM;
        totaltime(count,1) = KMS_output{mm}.totaltime;
        count = count+1;
    end
end

% Percent covered (out of Nmc).
cov = find(KMS_CI_feas(:,1) <= theta_true(component) &  KMS_CI_feas(:,2) >= theta_true(component));
percent_cov = size(cov,1)/Nmc;

% Coverage of lower/upper bound
covL =  find(KMS_CI_feas(:,1) <= Identification_region(component,1) &  KMS_CI_feas(:,2) >= Identification_region(component,1));
covU =  find(KMS_CI_feas(:,1) <= Identification_region(component,2) &  KMS_CI_feas(:,2) >= Identification_region(component,2));
percent_covL = size(covL,1)/Nmc;
percent_covU = size(covU,1)/Nmc;


% Report
med_KMS
count_flagfeas
count_flagLB
count_flagUB
mean_c = mean(c_theta)
mean_tEAM = mean(time_EAM)
mean_tt = mean(totaltime)
percent_cov
percent_covL
percent_covU












