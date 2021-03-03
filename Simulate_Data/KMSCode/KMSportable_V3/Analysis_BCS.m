% This code is the analog for Analysis.m for the BCS code.

% Vectorize the confidence itnerval
BCS_CI = zeros(Nmc,2);
for mm = 1:Nmc
    BCS_CI(mm,:) = BCS_confidence_interval{mm};
end

% Find runs that failed to converge and count
count_fail = 0;
ind_fail = zeros(Nmc,1);
total_time = zeros(Nmc,1);
for mm = 1:Nmc
    if BCS_output{mm}.flag_conv == 0
        count_fail = count_fail+1;
        ind_fail(mm,1) = 1;
    else
        total_time(mm,1) = BCS_output{mm}.totaltime;
    end
end
total_time(ind_fail == 1) = [];

% BCS results
med_BCS = median(BCS_CI(ind_fail==0,:));

% oUTPUT
count_fail
med_BCS
avg_minutes = mean(total_time)/(60)