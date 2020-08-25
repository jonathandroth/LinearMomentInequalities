%Calculate Linear IV regressions reported in Table V of Angrist Kreueger
%1991.  Here Consider same specifications as in Staiger and Stock (1997)
clear
rng(1);

%Specify whether to use laptop or server sepcifications
laptop=0;
if laptop==1
    chain_draws=10^3;
    %parpool('local',2)
else
    chain_draws=10^6;
    parpool('local',24)
end
size_sim_draws=10^5;
tau_grid_count=500;
D_draws=10^4;
coverage=.95;
%Specify whether to calculate plug-in confidence set
PI_calc=1;
%specify number of points to use in each dimension of beta for calculation
%of plug-in confidence set
PI_grid_points=5;
endog_reg=2;

beta_lower_bound=-.5;
beta_upper_bound=.5;

count=1;
while count<endog_reg
    beta_lower_bound(end+1,1)=-2;
    beta_upper_bound(end+1,1)=2;
    count=count+1;
end



%Use homoskedastic or heteroskedastic version of robust statistics
heterosked_robust=1;

%Choose which speficiation to use
spec=2;

load AK_3039_analysis_dataset;
st_dataset(end,:)=[];

AGEQ=double(st_dataset(1:end,2));
AGEQSQ=double(st_dataset(1:end,28));
EDUC=double(st_dataset(1:end,4));
ENOCENT=double(st_dataset(1:end,5));
ESOCENT=double(st_dataset(1:end,6));
LWKLYWGE=double(st_dataset(1:end,8));
MARRIED=double(st_dataset(1:end,9));
MIDATL=double(st_dataset(1:end,10));
MT=double(st_dataset(1:end,11));
NEWENG=double(st_dataset(1:end,12));
RACE=double(st_dataset(1:end,18));
SMSA=double(st_dataset(1:end,19));
SOATL=double(st_dataset(1:end,20));
WNOCENT=double(st_dataset(1:end,23));
WSOCENT=double(st_dataset(1:end,24));
YR20_28=double(st_dataset(1:end,29:37));
QTR120_129=double(st_dataset(1:end,43:52));
QTR220_229=double(st_dataset(1:end,53:62));
QTR320_329=double(st_dataset(1:end,63:72));
QTR1_3=double(st_dataset(1:end,39:41));
v17=double(st_dataset(1:end,16));

SOB_index=[1 2 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 ...
    32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55]; %drop 56 to avoid multicollinartity
SOB_dummies=zeros(size(QTR1_3,1),length(SOB_index));
for n=1:length(SOB_index)
    SOB_dummies(:,n)=double(v17==SOB_index(n));
end

JR_HS=(EDUC>=8);
Sophmore=(EDUC>=10);
Junior=(EDUC>=11);
HS=(EDUC>=12);
JR_College=(EDUC>=14);
College=(EDUC>=16);
Masters=(EDUC>=18);



%%

%Outcome variable is wages, endogenous variable is education, Z is
%instruments, and W are controls
Y=LWKLYWGE;
%X=[EDUC JR_HS];
%X=[EDUC Junior];
X=[EDUC HS];
%X=[EDUC JR_College];
%X=[EDUC HS JR_College];
%X=[EDUC HS College];
%X=[EDUC EDUC.^2];
%X=[EDUC  EDUC.*(EDUC>=12)];
%X=[EDUC HS EDUC.*(EDUC>=12)];
if endog_reg==1
    X=EDUC;
elseif endog_reg==2
    X=[EDUC HS];
elseif endog_reg==3
    X=[EDUC HS JR_College];
elseif endog_reg==4
    X=[EDUC HS JR_College College];
elseif endog_reg==5
    X=[EDUC HS JR_College College Masters];
end


%Collect exogenous control variables and instruments.  Have checked that
%results in 1-3 match corresponding results in Staiger and Stock.
if spec==1
    Z=QTR1_3;
    W=[ones(length(Y),1) RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT YR20_28];
elseif spec==2
    W=[ones(length(Y),1) RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT YR20_28];
    Z=[QTR120_129 QTR220_229 QTR320_329];
elseif spec==3
    W=[ones(length(Y),1) RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT AGEQ AGEQSQ YR20_28];
    Z=[QTR120_129 QTR220_229 QTR320_329];
    Z(:,end)=[];
    Z(:,end-1)=[];
else
    W=[ones(length(Y),1) RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT AGEQ AGEQSQ YR20_28 SOB_dummies];
    Z=[QTR120_129 QTR220_229 QTR320_329];
    Z(:,end)=[];
    Z(:,end-1)=[];
    Z=[Z repmat(SOB_dummies,1,3).*repmat(QTR1_3,1,size(SOB_dummies,2))];
end

%Drop rows with missing observations
data=[Y X Z W];
nanrows=isnan(sum(data,2));
Y=Y(nanrows==0,:);
X=X(nanrows==0,:);
Z=Z(nanrows==0,:);
W=W(nanrows==0,:);

clear data SOB_index AGEQ AGEQSQ EDUC ENOCENT ESOCENT LWKLYWGE MARRIED QTR220_229 QTR1_3 st_dataset
clear MIDATL MT NEWENG RACE SMSA SOATL WNOCENT WSOCENT YR20_28 QTR120_129 QTR320_329 SOB_dummies v17
clear College JR_College JR_HS Junior Sophmore Masters nanrows

%Run as 2SLS with included exogenous variables
X_tilde=[X W];
Z_tilde=[Z W];

ZZ=Z_tilde'*Z_tilde;
ZZ_inv=ZZ^-1;
ZX=Z_tilde'*X_tilde;
XZ=X_tilde'*Z_tilde;
ZY=Z_tilde'*Y;

beta_tilde=(XZ*ZZ_inv*ZX)^-1*XZ*ZZ_inv*ZY; %Matches estimate from AK code and paper

WW=W'*W;
WZ=W'*Z;
WX=W'*X;
WY=W'*Y;

Y=Y-W*WW^-1*WY;
X=X-W*WW^-1*WX;
Z=Z-W*WW^-1*WZ;

ZZ=Z'*Z;
ZZ_inv=ZZ^-1;
ZX=Z'*X;
XZ=X'*Z;
ZY=Z'*Y;

instruments=size(Z,2);
obs=length(Y);
Omega=NaN((size(X,2)+1)*instruments);
hmat=repmat(Y,1,instruments).*Z;
Omega(1:instruments,1:instruments)=cov(hmat);
clear W hmat X_tilde Z_tilde
for m1=0:size(X,2)
    if m1==0
        hmat1=repmat(Y,1,instruments).*Z;
    else
        hmat1=repmat(X(:,m1),1,instruments).*Z;
    end
    for m2=1:size(X,2)
        hmat2=repmat(X(:,m2),1,instruments).*Z;
        Omega(m1*instruments+1:(m1+1)*instruments,m2*instruments+1:(m2+1)*instruments)...
            =1/obs*hmat1'*hmat2-mean(hmat1)*mean(hmat2)';
        Omega(m2*instruments+1:(m2+1)*instruments,m1*instruments+1:(m1+1)*instruments)...
            =1/obs*(hmat1'*hmat2-mean(hmat1)*mean(hmat2)')';
    end
end

clear hmat1 hmat2 X Y Z;
Sigma_gg=@(b) kron([1;-b],eye(instruments))'*Omega*kron([1;-b],eye(instruments));
Sigma_tg=@(b) Omega(instruments+1:end,1:instruments)...
    -Omega(instruments+1:end,instruments+1:end)*kron(b,eye(instruments));
ZX_vec=reshape(ZX,length(ZX(1:end)),1);
g=@(b) ZY-ZX*b;
D=@(b) reshape(ZX_vec-Sigma_tg(b)*Sigma_gg(b)^-1*g(b),instruments,size(ZX,2));
S=@(b) 1/obs*g(b)'*Sigma_gg(b)^-1*g(b);

K_draws=sum(randn(size_sim_draws,size(ZX,2)).^2,2);
J_draws=sum(randn(size_sim_draws,instruments-size(ZX,2)).^2,2);
rk_grid=[0:.1:10 11:1:100 110:10:1000].^2;

QCLR_crit=zeros(length(rk_grid),1);
for m=1:length(rk_grid)
    r=rk_grid(m);
    QCLR_crit(m,1)=quantile(.5*(K_draws+J_draws-r+...
        ((K_draws+J_draws+r).^2-4*J_draws*r).^.5),coverage);
end

chain_store=zeros(chain_draws,size(ZX,2));
S_store=ones(chain_draws,1)*1000;
K_store=ones(chain_draws,1)*1000;
rk_store=zeros(chain_draws,1);
QCLR_stat_store=zeros(chain_draws,1);
QCLR_crit_store=zeros(chain_draws,1);
parfor m=2:chain_draws
    if round(m*100/chain_draws)==m*100/chain_draws
        m/chain_draws
    end
    
    proposal=(rand(size(ZX,2),1)-.5).*(beta_upper_bound-beta_lower_bound)+(beta_upper_bound+beta_lower_bound)/2;
    
    
    A=kron([1;-proposal],eye(instruments));
    Sigma_temp=A'*Omega*A;
    Sigma_temp_inv=Sigma_temp^-1;
    Sigma_tg_temp=Omega(instruments+1:end,1:instruments)...
        -Omega(instruments+1:end,instruments+1:end)*kron(proposal,eye(instruments));
    g_temp=g(proposal);
    S_val=1/obs*g_temp'*Sigma_temp_inv*g_temp;
    S_store(m,1)=S_val;
    D_temp=reshape(ZX_vec-Sigma_tg_temp*Sigma_temp_inv*g_temp,instruments,size(ZX,2));
    
    K_store(m,1)=1/obs*g_temp'*Sigma_temp_inv*D_temp*(D_temp'*Sigma_temp_inv*D_temp)...
        ^-1*D_temp'*Sigma_temp_inv*g_temp;
    chain_store(m,:)=proposal;
    J_temp=S_store(m,1)-K_store(m,1);
    
    if size(ZX,2)==1;
        D_temp=D_temp';
    end
    Sigma_D_temp=Omega(instruments+1:end,instruments+1:end)-Sigma_tg_temp...
        *Sigma_temp_inv*Sigma_tg_temp';
    Sigma_D_temp_inv_rt=Sigma_D_temp^-.5;
    D_normalized=reshape(Sigma_D_temp_inv_rt*D_temp(1:end)',instruments,size(ZX,2));
    rk_store(m,1)=1/obs*min(eig(D_normalized'*D_normalized));
    QCLR_stat_store(m,1)=.5*(K_store(m,1)+J_temp-rk_store(m,1)+...
        ((K_store(m,1)+J_temp+rk_store(m,1)).^2-4*J_temp*rk_store(m,1)).^.5);
    QCLR_crit_store(m,1)=interp1q(rk_grid,QCLR_crit,rk_store(m,1));
end

chain_truncate_S=chain_store(S_store<chi2inv(coverage,instruments),:);
chain_truncate_K=chain_store(K_store<chi2inv(coverage,size(ZX,2)),:);
chain_truncate_QCLR=chain_store(QCLR_stat_store<QCLR_crit_store,:);
chain_truncate_JK=chain_store((K_store<chi2inv(1-(1-coverage)*.8,size(ZX,2)))&((S_store-K_store)<chi2inv(1-(1-coverage)*.2,instruments-size(ZX,2))),:);

a_grid=[0:.01:1];
identity_vec=zeros(chain_draws,1);
for n=1:length(a_grid);
    cutoff(n)=quantile((1-a_grid(n))*K_draws+a_grid(n)*J_draws,coverage);
    chain_truncate_LC=chain_store((1-a_grid(n))*K_store+a_grid(n)*S_store<cutoff(n),:);
    LC_vol(n,1)=length(chain_truncate_LC(:,1))/chain_draws;
    identity_vec((1-a_grid(n))*K_store+a_grid(n)*S_store<cutoff(n))=1;
    for m=1:size(ZX,2)
        LC_CI(n,:,m)=[min(chain_truncate_LC(:,m)) max(chain_truncate_LC(:,m))];
    end
end
[C,I] = min(LC_vol);
min_vol_a=a_grid(I);
min_a_cutoff=quantile((1-min_vol_a)*K_draws+min_vol_a*J_draws,coverage);
chain_truncate_min_a_LC=chain_store((1-min_vol_a)*K_store+min_vol_a*S_store<min_a_cutoff,:);


if PI_calc==1
    
    vecmat=beta_lower_bound(1):((beta_upper_bound(1)-beta_lower_bound(1))/(PI_grid_points-1)):beta_upper_bound(1);
    
    for n=2:length(beta_lower_bound)
        vecmat(n,:)=beta_lower_bound(n):((beta_upper_bound(n)-beta_lower_bound(n))/(PI_grid_points-1)):beta_upper_bound(n);
    end
    
    if endog_reg==1
        M1=vecmat;
        bmat=M1(1:end)';
        load Powerblock_dk1_dj29_1Msims
    elseif endog_reg==2
        [M1, M2]=ndgrid(vecmat(1,:),vecmat(2,:));
        bmat=[M1(1:end)' M2(1:end)'];
        load Powerblock_dk2_dj28_1Msims
    elseif endog_reg==3
        [M1, M2, M3]=ndgrid(vecmat(1,:),vecmat(2,:),vecmat(3,:));
        bmat=[M1(1:end)' M2(1:end)' M3(1:end)'];
        load Powerblock_dk3_dj27_1Msims
    elseif endog_reg==4
        [M1, M2, M3, M4]=ndgrid(vecmat(1,:),vecmat(2,:),vecmat(3,:),vecmat(4,:));
        bmat=[M1(1:end)' M2(1:end)' M3(1:end)' M4(1:end)'];
        load Powerblock_dk4_dj26_1Msims
    elseif endog_reg==5
        [M1, M2, M3, M4, M5]=ndgrid(vecmat(1,:),vecmat(2,:),vecmat(3,:),vecmat(4,:),vecmat(5,:));
        bmat=[M1(1:end)' M2(1:end)' M3(1:end)' M4(1:end)' M5(1:end)'];
        load Powerblock_dk5_dj25_1Msims
    end
    clear M1 M2 M3 M4 M5
    
    
    chain_truncate_all_LC=chain_store(identity_vec==1,:);
    D_noise2=randn(D_draws,instruments*size(ZX,2));
    PIRMS_test=zeros(size(chain_truncate_all_LC,1),1);
    location_store=zeros(size(chain_truncate_all_LC,1),1);
    a_store=zeros(size(chain_truncate_all_LC,1),1);
    tic;
    parfor m=1:size(chain_truncate_all_LC,1)
    %for m=1:size(chain_truncate_all_LC,1)
        location_store(m)=1;
        m/size(chain_truncate_all_LC,1)
        beta0=chain_truncate_all_LC(m,:)';
        A=kron([1;-beta0],eye(instruments));
        Sigma_temp=A'*Omega*A;
        Sigma_temp_inv=Sigma_temp^-1;
        Sigma_temp_inv_rt=real(Sigma_temp_inv^.5);
        Sigma_tg_temp=Omega(instruments+1:end,1:instruments)...
            -Omega(instruments+1:end,instruments+1:end)*kron(beta0,eye(instruments));
        Sigma_D_temp=Omega(instruments+1:end,instruments+1:end)-Sigma_tg_temp...
            *Sigma_temp_inv*Sigma_tg_temp';
        Sigma_D_temp_rt=real(Sigma_D_temp^.5);
        g_temp=g(beta0);
        S_temp=1/obs*g_temp'*Sigma_temp_inv*g_temp;
        D_temp=reshape(ZX_vec-Sigma_tg_temp*Sigma_temp_inv*g_temp,instruments,size(ZX,2));
        K_temp=1/obs*g_temp'*Sigma_temp_inv*D_temp*(D_temp'*Sigma_temp_inv*D_temp)...
            ^-1*D_temp'*Sigma_temp_inv*g_temp;
        J_temp=S_temp-K_temp;
        location_store(m)=2;
        %Just using plug-in estimate: can likely do something more refined
        %mu_hat=obs^-.5*D_temp;
        
        %Use natural generalization of positive-part estimator
        V=zeros(instruments,instruments);
        for n=1:size(D_temp,2)
            V=V+Sigma_D_temp((n-1)*instruments+1:n*instruments,(n-1)*instruments+1:n*instruments);
        end
        if max(max(isnan(V)))==0&&max(max(isinf(V)))==0
            if rank(V)==size(V,2)
                V=(1/size(D_temp,2)*V)^-1;
            else
                V=eye(size(V,2));
            end
        else
            V=eye(size(V,2));
        end
        R=obs^-1*D_temp'*V*D_temp;
        location_store(m)=3;
        Bias_mat=zeros(size(D_temp,2));
        for n1=1:size(D_temp,2)
            for n2=1:size(D_temp,2)
                Bias_mat(n1,n2)=trace(V*Sigma_D_temp((n1-1)*instruments+1:n1*instruments,(n2-1)*instruments+1:n2*instruments));
            end
        end
        R_U=R-Bias_mat;
        [V L]=eig(R_U);
        L((L<0)&(isnan(L)==0))=0;
        R_U=V*L*V';
        location_store(m)=4;
        if max(max(isnan(R)))==0&&max(max(isinf(R)))==0
            if rank(R)==size(R,2)
                mu_D_hat=obs^-.5*D_temp*real(R^-.5*R_U^.5);
            else
                mu_D_hat=obs^-.5*D_temp;
            end
        else
            mu_D_hat=obs^-.5*D_temp;
        end
        mu_D_hat=real(mu_D_hat);
        location_store(m)=5;
        D_noise=D_noise2*Sigma_D_temp_rt;
        tau_J=zeros(D_draws,size(bmat,1));
        tau_K=zeros(D_draws,size(bmat,1));
        
        %Calculate m under different alternatives
        m_hat=zeros(instruments,size(bmat,1));
        for n=1:size(bmat,1)
            b_diff=bmat(n,:)'-beta0;
            vec_mu_hat=(eye(size(ZX,2)*instruments)-kron(b_diff',Sigma_tg_temp*Sigma_temp_inv))\mu_D_hat(1:end)';
            mu_hat=reshape(vec_mu_hat,instruments,size(ZX,2));
            m_hat(:,n)=mu_hat*b_diff;
        end
        
        %Calculate tau_J and tau_K for simulated draws of D
        for n=1:D_draws
            D_base=reshape(D_noise(n,:),instruments,size(ZX,2));
            D_draw=D_base+mu_D_hat;
            temp=D_draw'*Sigma_temp_inv*D_draw;
            if max(max(isnan(temp)))==0&&max(max(isinf(temp)))==0
                if rank(D_draw'*Sigma_temp_inv*D_draw)==size(D_draw,2);
                    D_norm=real((D_draw'*Sigma_temp_inv*D_draw)^-.5*D_draw'*Sigma_temp_inv);
                else
                    D_norm=D_draw'*0;
                end
            else
                D_norm=D_draw'*0;
            end
            tau_K(n,:)=sum((D_norm*m_hat).^2,1);
            tau_J(n,:)=sum((Sigma_temp_inv_rt*m_hat).^2,1)-tau_K(n,:);
        end
        location_store(m)=6;
        cell_mean_TK=zeros(size(bmat,1),tau_grid_count+1);
        cell_mean_TJ=zeros(size(bmat,1),tau_grid_count+1);
        hist_est=zeros(size(bmat,1),tau_grid_count+1);
        for n=1:size(bmat,1)
            hist_est(n,:)=ones(1,tau_grid_count+1)/(tau_grid_count+1);
            tmp=[min(tau_K(:,n)) quantile(tau_K(:,n),tau_grid_count) max(tau_K(:,n))];
            cell_mean_TK(n,:)=(tmp(1:end-1)+tmp(2:end))/2;
            cell_mean_TJ(n,:)=sum((Sigma_temp_inv_rt*m_hat(:,n)).^2,1)-cell_mean_TK(n,:);
        end
        location_store(m)=7;
        
        powertemp=zeros(size(bmat,1),size(a_grid,2));
        for na=1:size(a_grid,2)
            cond_power=interp2(TJ,TK,powerblock(:,:,na),cell_mean_TJ,cell_mean_TK,'linear',1);
            avg_power=sum(cond_power.*hist_est,2);
            powertemp(:,na)=avg_power;
        end
        location_store(m)=8;
        power_max=max(powertemp,[],2);
        max_deficiency=max(repmat(power_max,1,size(powertemp,2))-powertemp,[],1);
        a=a_grid(max_deficiency<=min(max_deficiency)+10^-5);
        
        if isempty(a)==0
            a_PIRMS=a(end);
        else
            a_PIRMS=1;
        end
        PIRMS=K_temp+a_PIRMS*J_temp;
        PIRMS_test(m)=PIRMS>interp1(a_grid,a_crit,a_PIRMS);
        a_store(m)=a_PIRMS;
        location_store(m)=9;
    end
    chain_truncate_PI=chain_truncate_all_LC(PIRMS_test==0,:);
end
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers
end
total_PI_runtime=toc*poolsize;
average_PI_runtimes=toc/size(chain_truncate_all_LC,1)*poolsize;


clear bmat D_noise2 K_draws J_draws TJ TK
if laptop==1
    delete(gcp)
else
    if endog_reg==1
        save AK_CLC_1X_5points_1M_draws
    elseif endog_reg==2
        save AK_CLC_2X_5points_1M_draws
    elseif endog_reg==3
        save AK_CLC_3X_5points_1M_draws
    elseif endog_reg==4
        save AK_CLC_4X_5points_1M_draws
    elseif endog_reg==5
        save AK_CLC_5X_5points_1M_draws
    end
    
    exit
end