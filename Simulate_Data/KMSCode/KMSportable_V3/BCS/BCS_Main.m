function [BCS_CI,BCS_output] = BCS_Main(W,theta_true,component,selp,n,B,seed,alpha)

%% Preliminary

coordinate = component;

%% Set parameters for simulation
%alpha = 0.1; % significance level
%display_info = 1; % verbose on/off

%% Detect Simulation Type
%directory_name = pwd; % record directory name
% directory should be named SOMENAME_x, where x = denotes coordinate of interest for inference (1 or 2).
%coordinate = str2double(directory_name(end)); % Grab x, last character from directory name.
%assert(coordinate==1 || coordinate==2);

%% Primitives of DGP
%theta0 = [0.3,0.5,0,0,0.05]; % true parameter of interest. In the paper is theta0 = (theta1,theta2,beta2,beta3,beta1).
theta0 = theta_true.';
dimSX = size(theta0,2)-1; % dimension of X.
PX = ones(dimSX,1)/dimSX; % P(X=x), X is discrete uniform.
prob10 = selp; %prob10 = 0.6; % prob of P(A_1=1,A_2=0) when there is multiplicity.

%% Grid of values of theta according to H0;
theta_delta = 0.005;
if coordinate == 1
    theta_H0_list = (0.005:theta_delta:0.995)';
else %coordinate = 2.
    theta_H0_list = (0.005:theta_delta:0.995)';
end

%% Define boundaries of the Id Set. These numbers were numerically calculated by Ivan
% if coordinate == 1
%     trueH0  = (theta_H0_list>=0.23).*(theta_H0_list<=0.36);
% else % coordinate = 2.
%     trueH0  = (theta_H0_list>=0.455).*(theta_H0_list<=0.555);
% end

% sample and bootstrap sample sizes;
%n = 1000; % sample size.
R = B; %R = 300; % bootstrap sample size, i.e., the simulated quantiles are the result of R simulated values.
%MC_size = 2000; % MC size, i.e., empirical coverage is the result of MC_size tests.

%% Tuning parameters for inference;
bn_list =  round([n/4;n^(2/3)]); % subsampling tuning parameter.
C_list = 1; %C_list = [0.8;1]; % constant in the paper;
kappa_list = C_list*sqrt(log(n));  % GMS Thresholding parameter, preferred choice Eq. 4.4 in Andrews and Soares (2010).

%% Parameters for optimization;
tolerance = 0.0001; % define what we mean by "close enough".

% options for minimization.
options = optimset('Display','off','Algorithm','interior-point'); % (this is slower).
% options = optimset('Display','off','Algorithm','active-set'); % (this is faster).

%random_seed_start = 1; % Arbitrary start for the seed;
%random_seed_list  = random_seed_start+(1:1:MC_size); % list or random seeds.

%disp('Save all information in the BASE directory');
% Create directory to save this information;
%mkdir results;
%cd results
matrix_reject_test_DR = NaN(1,size(theta_H0_list,1),size(kappa_list,1));
matrix_reject_test_PR = NaN(1,size(theta_H0_list,1),size(kappa_list,1));
matrix_reject_test_MR = NaN(1,size(theta_H0_list,1),size(kappa_list,1));
%matrix_reject_test_SS = NaN(MC_size,size(theta_H0_list,1),size(bn_list,1));
%save results matrix_reject_test_*;
%cd ..

%mkdir BASE
%cd BASE
%save Data.mat alpha display_info coordinate theta0 dimSX PX prob10 theta_H0_list trueH0 n R bn_list C_list kappa_list tolerance options MC_size random_seed_list;
%cd ..


%% Main BCS Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Test program: This programs runs all tests for all MC simulations. It
%   should be ran after the preliminary program.
%   In practice, we paralellized the 2000 MC runs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load basic information;
%cd BASE
%load Data.mat MC_size display_info
%cd ..

% MC loop
%for MC_index = 1 : MC_size
%disp(['* MC index = ',num2str(MC_index),' out of ',num2str(MC_size)]);

% Load basic information;
%cd BASE
%load Data.mat theta_H0_list
%cd ..

% Loop over all values of theta of interest under H0

% Set seed (new one for each # bootstrap)
try

    %% Lower bound
    for theta_index = 1:size(theta_H0_list,1)
        stream=RandStream('mlfg6331_64','Seed',seed);
        RandStream.setGlobalStream(stream);
        stream.Substream = theta_index + B*10^3;
        % check: is simulation already done?

        % go into results directory
        %cd results;
        % only load one of the results matrices;
        %load results matrix_reject_test_MR;
        %cd ..

        %if isnan(matrix_reject_test_MR(MC_index,theta_index,1))==1
            % this happens if not done: does the simulation.
        %disp(['** theta_index = ',num2str(theta_index),' out of ',num2str(size(theta_H0_list,1))]);

        %tic; % start the time;

        %% Load basic information and initialize matrices
        %cd BASE
        %load Data.mat;
        %cd ..

        % fix the simulation seed depending on MC_index;
        %estado_random = random_seed_list(MC_index); % initialize random seed using the saved info
        %rand('state',estado_random); %#ok<RAND> % use the seed as intial state.

        % simulate basic random variables we use to generate data;
        p = 2*dimSX; % total inequalities;
        k = 3*dimSX; % total (in)equalties;

        %baseDatas = rand(n,4); % this information determines all random draws.

        % Random draws for BP/Projection test. We don't use it as we did not implement BP.
        W1 = norminv(rand(R,k));

        % Random draws for RC, RS, and MR tests
        W2_AA = norminv(rand(R,n));

        % Random draws for SS test
    %     W3_SS = NaN(R,max(bn_list));
    %     for r=1:R
    %         aux = randperm(n); % randomly permute the n observations;
    %         W3_SS(r,:) = aux(1:max(bn_list)); % choose the first bn observations;
    %     end
    %     clear aux*;

        [~,chi] = S_function([],[]); % Get degree of homogeneity from definition of S1;

        %% Initialize matrices

        % Initialize matrices for SS
        %minQn_SS = NaN(size(bn_list,1),size(W3_SS,1));
        cn_SS = NaN(1,size(bn_list,1));

        % Initialize matrices for GMS: BP, RC, RS
        minQn_DR = NaN(size(kappa_list,1),size(W2_AA,1));
        minQn_PR = NaN(size(kappa_list,1),size(W2_AA,1));
        minQn_MR = NaN(size(kappa_list,1),size(W2_AA,1));
        cn_DR  = NaN(1,size(kappa_list,1));
        cn_PR  = NaN(1,size(kappa_list,1));
        cn_MR  = NaN(1,size(kappa_list,1));

        %% Simulate data
        %epsilons = baseDatas(:,1:2); % epsilon in the model in section 5;
        %multiple = baseDatas(:,3); % determines how multiplicity is resolved;

        % X denotes the market type indicator;
        %X = zeros(size(baseDatas(:,4)));
        %for j=1:dimSX
        %    X = X + (j-1)*(baseDatas(:,4)>=(j-1)/dimSX).*(baseDatas(:,4)<j/dimSX) ;
        %end
        %betas0 = [0,theta0(3:end)]; % vector indicates beta_q for q=1,...,d_X

        % Initialize matrices that will contain both entry decision and market type
    %     dataP11 = NaN(size(baseDatas,1),size(PX,1));
    %     dataP10 = NaN(size(baseDatas,1),size(PX,1));
    %     for j=1:dimSX
    % 
    %         % Entry decision that indices {A_1=1,A_2=1}
    %         Entry11_aux = (epsilons(:,1) > theta0(1)-betas0(j)).*(epsilons(:,2) >  theta0(2)-betas0(j)) ;
    %         % should have mean = (1-(theta0(1)-betas0(j)))*(1-(theta0(2)-betas0(j)));
    % 
    %         % Entry decision that indices {A_1=1,A_2=0}
    %         Entry10_aux = (epsilons(:,1) > theta0(1)-betas0(j)).*(epsilons(:,2) <= theta0(2)-betas0(j)) + (epsilons(:,1) <= theta0(1)-betas0(j)).*(epsilons(:,2) <= theta0(2)-betas0(j)).*(multiple>prob10);
    %         % should have mean = (theta0(1)-betas0(j))*(theta0(2)-betas0(j))*(1-prob10) + (1-(theta0(1)-betas0(j)))*(theta0(2)-betas0(j));
    % 
    %         % entry (positive or zero) and market type (integer) information (It assumes P(X=x) is known)
    %         dataP11(:,j) = Entry11_aux.*(X == j-1)/PX(j); % indicates {A_1=1,A_2=1} and {X=j-1}
    %         dataP10(:,j) = Entry10_aux.*(X == j-1)/PX(j); % indicates {A_1=1,A_2=0} and {X=j-1}
    %     end

        % collect entry data in a single matrix;
        data =  W; %data = [dataP11,dataP10];
        dataP11 = W(:,1:dimSX);
        dataP10 = W(:,dimSX+1:end);

        % determine the value of theta of interest according to H0;
        theta_H0 = theta_H0_list(theta_index);

        % determines the value of the other thetas;
        if coordinate==1
            theta_other0 = theta0(2:end);
        else % coordinate = 2
            theta_other0 = [theta0(1),theta0(3:end)];
        end

        % indicates if H0 is true or not;
        %theta_in_IS = trueH0(theta_index);

        %% Compute test statistic for all tests, profile test statistic.

        % Choose starting values for minimization;
        starting_values = NaN(2,size(theta_other0,2)); % initialize;

        % first starting value: true parameter value (unfeasible in practice)
        starting_values(1,:) = theta_other0';

        % second starting value: good guess based on the theory (feasible in practice)
        starting_values(2,1) = 1 - mean(dataP11(:,1))/(1-theta_H0);
        % the remaining starting values are randomly chosen;
        for a=2:size(starting_values,2)
            if theta_other0(a)>0
                starting_values(2,a) = -(1/2)*(2-theta_H0-starting_values(2,1))+(1/2)*sqrt( (starting_values(2,1) - theta_H0)^2 + 4*mean(dataP11(:,a)) );
            else
                starting_values(2,a)=0;
            end
        end

        % Pick a large value as initial function value for the minimization
        min_value = 10^10;

        % for minimization we need the following constrains: Aineq*theta0?Bineq and (0,...,0)<=theta0<=(1,...1)
        Aineq = zeros(2*(size(theta0,2)-2),size(theta0,2)-1);
        Aineq(1:size(theta0,2)-2,2:size(theta0,2)-1) = eye(size(Aineq(1:size(theta0,2)-2,2:size(theta0,2)-1)));
        Aineq(size(theta0,2)-1:end,2:size(theta0,2)-1) = eye(size(Aineq(size(theta0,2)-1:end,2:size(theta0,2)-1)));
        Aineq(size(theta0,2)-1:end,1) = -1*ones(size(Aineq(size(theta0,2)-1:end,1)));
        Bineq = zeros(2*(size(theta0,2)-2),1);
        Bineq(1:(size(theta0,2)-2)) = theta_H0 * ones(size(Bineq(1:(size(theta0,2)-2))));

        min_outcomes = []; % This matrix will collect results of minimization;

        % solve numerical minimization for all starting values;
        for s=1:size(starting_values,1)
            [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_function(x,theta_H0,coordinate,data,p),starting_values(s,:),Aineq,Bineq,[],[],zeros(size(theta_other0)),ones(size(theta_other0)),[],options);
            % check whether minimization is successful and reduced value
            if Qn_aux  < min_value && bandera >= 1
                Qn_minimizer = theta_aux;
                min_value = Qn_aux;
            end
            % if minimization is successful, collect minimizer and its value;
            if bandera >= 1
                min_outcomes = [min_outcomes;[min_value,Qn_minimizer]]; %#ok<AGROW>
            end
        end
        % at the end of the minimizations, we should have a minimizer
        minQn = min_value;

        % Collect minimizers to estimate the indetified set (used for DR)
        min_outcomes = uniquetol2(min_outcomes,tolerance,'rows'); % set of minimizers;
        In_identified_set_hat = min_outcomes(min_outcomes(:,1)<=min(min_outcomes(:,1)) + tolerance , 2:end); % estimator of id set

        %% BP, i.e., projection.
        % By definition, this method requires gridsearch. It is feasible in simpler problems but not in this model and using Matlab.

        %% Subsampling test
    %     for bn_type=1:size(bn_list,1)
    % 
    %         for r=1:R
    %             % collect relevant data for SS minimization;
    %             data_SS = data(W3_SS(r,1:bn_list(bn_type)),:);
    % 
    %             % Starting values for SS minimization;
    %             starting_values =  uniquetol2([starting_values;Qn_minimizer],tolerance,'rows');
    % 
    %             % Pick a large value as initial function value for the minimization;
    %             min_value_SS = 10^10;
    % 
    %             % solve numerical minimization for all starting values;
    %             for s=1:size(starting_values,1)
    %                 [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_function(x,theta_H0,coordinate,data_SS,p),starting_values(s,:),Aineq,Bineq,[],[],zeros(size(theta_other0)),ones(size(theta_other0)),[],options);
    %                 % check whether minimization is successful and reduced value
    %                 if Qn_aux  < min_value_SS && bandera >= 1
    %                     minimizer = theta_aux;
    %                     min_value_SS = Qn_aux;
    %                 end
    %             end
    %             minQn_SS(bn_type,r) = min_value_SS;
    %         end
    % 
    %         % Compute SS critical value
    %         cn_SS(bn_type) = quantile(minQn_SS(bn_type,:),1-alpha);
    %     end
    %     % SS Test: Accept decision;
    %     accept_test_SS = repmat(minQn,1,size(bn_list,1)) <= cn_SS;

        % Starting values for PR minimization;
        starting_values =  uniquetol2([starting_values;Qn_minimizer],tolerance,'rows');

        %% MR test (and also DR and PR)
        for kappa_type =1:size(kappa_list,1)
            parfor r=1:R

                %- Step 1: Discard resampling method

                % Pick a large value as initial function value for the minimization;
                min_value_DR = 10^10;

                % Minimize but restricted to points of the estimated indentified set;
                for IDsetHat_index = 1:size(In_identified_set_hat,1)
                    minQn_DR_aux = min(min_value_DR,Qn_MR_function(In_identified_set_hat(IDsetHat_index,:),theta_H0,coordinate,data,kappa_list(kappa_type),p,k,W2_AA(r,:),1));
                end

                % compute simulated DR criterion function
                minQn_DR(kappa_type,r) = minQn_DR_aux;

                %- Step 2: Penalized resampling method;

                % Pick a large value as initial function value for the minimization;
                min_value_PR = 10^10;
                % check whether minimization is successful and reduced value
                for s=1:size(starting_values,1)
                    [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_MR_function(x,theta_H0,coordinate,data,kappa_list(kappa_type),p,k,W2_AA(r,:),2),starting_values(s,:),Aineq,Bineq,[],[],zeros(size(theta_other0)),ones(size(theta_other0)),[],options);
                    if Qn_aux  < min_value_PR && bandera >= 1
                        minimizer = theta_aux;
                        min_value_PR = Qn_aux;
                    end
                end

                % compute simulated PR criterion function
                minQn_PR(kappa_type,r) = min_value_PR ;

                %- Step 3: combine PR and DR to get MR
                minQn_MR(kappa_type,r) = min(minQn_DR(kappa_type,r),minQn_PR(kappa_type,r));
            end
            % Compute critical values
            cn_DR(kappa_type) = quantile(minQn_DR(kappa_type,:),1-alpha);
            cn_PR(kappa_type) = quantile(minQn_PR(kappa_type,:),1-alpha);
            cn_MR(kappa_type) = quantile(minQn_MR(kappa_type,:),1-alpha);
        end

        % DR, PR, and MR Tests: Accept decisions
        accept_test_DR = repmat(minQn,1,size(kappa_list,1)) <= cn_DR;
        accept_test_PR = repmat(minQn,1,size(kappa_list,1)) <= cn_PR;
        accept_test_MR = repmat(minQn,1,size(kappa_list,1)) <= cn_MR;

        %% Finish simulation: save and display time;
        %tElapsed = toc;

        % Save newly found results
        %cd results
        % load the corresponding matrices;
        %load results
        % replace the NaN with the corresponding results;
        %matrix_reject_test_SS(MC_index,theta_index,:) = ones(size(accept_test_SS)) - accept_test_SS; %#ok<SAGROW>
        matrix_reject_test_DR(1,theta_index) = ones(size(accept_test_DR)) - accept_test_DR; %#ok<SAGROW>
        matrix_reject_test_PR(1,theta_index) = ones(size(accept_test_PR)) - accept_test_PR; %#ok<SAGROW>
        matrix_reject_test_MR(1,theta_index) = ones(size(accept_test_MR)) - accept_test_MR;
        %save results matrix_reject_test*
        %cd ..

        % Display overall time of the MCrun in minutes
        %if display_info == 1
        %    disp(['- MCrun done in ' num2str(tElapsed/60) ' mins.']);
        %end
        %fprintf('iter %9.4f of %9.4f done \n', [theta_index,size(theta_H0_list,1)])


        % If accept MR test, then we have found the lower bound.
        if matrix_reject_test_MR(1,theta_index) == 0
            break;
        end
    end

%% Upper bound
    for theta_index = size(theta_H0_list,1):-1:1
        stream=RandStream('mlfg6331_64','Seed',seed);
        RandStream.setGlobalStream(stream);
        stream.Substream = theta_index + B*10^3;
        % check: is simulation already done?

        % go into results directory
        %cd results;
        % only load one of the results matrices;
        %load results matrix_reject_test_MR;
        %cd ..

        %if isnan(matrix_reject_test_MR(MC_index,theta_index,1))==1
            % this happens if not done: does the simulation.
        %disp(['** theta_index = ',num2str(theta_index),' out of ',num2str(size(theta_H0_list,1))]);

        %tic; % start the time;

        %% Load basic information and initialize matrices
        %cd BASE
        %load Data.mat;
        %cd ..

        % fix the simulation seed depending on MC_index;
        %estado_random = random_seed_list(MC_index); % initialize random seed using the saved info
        %rand('state',estado_random); %#ok<RAND> % use the seed as intial state.

        % simulate basic random variables we use to generate data;
        p = 2*dimSX; % total inequalities;
        k = 3*dimSX; % total (in)equalties;

        %baseDatas = rand(n,4); % this information determines all random draws.

        % Random draws for BP/Projection test. We don't use it as we did not implement BP.
        W1 = norminv(rand(R,k));

        % Random draws for RC, RS, and MR tests
        W2_AA = norminv(rand(R,n));

        % Random draws for SS test
    %     W3_SS = NaN(R,max(bn_list));
    %     for r=1:R
    %         aux = randperm(n); % randomly permute the n observations;
    %         W3_SS(r,:) = aux(1:max(bn_list)); % choose the first bn observations;
    %     end
    %     clear aux*;

        [~,chi] = S_function([],[]); % Get degree of homogeneity from definition of S1;

        %% Initialize matrices

        % Initialize matrices for SS
        %minQn_SS = NaN(size(bn_list,1),size(W3_SS,1));
        cn_SS = NaN(1,size(bn_list,1));

        % Initialize matrices for GMS: BP, RC, RS
        minQn_DR = NaN(size(kappa_list,1),size(W2_AA,1));
        minQn_PR = NaN(size(kappa_list,1),size(W2_AA,1));
        minQn_MR = NaN(size(kappa_list,1),size(W2_AA,1));
        cn_DR  = NaN(1,size(kappa_list,1));
        cn_PR  = NaN(1,size(kappa_list,1));
        cn_MR  = NaN(1,size(kappa_list,1));

        %% Simulate data
        %epsilons = baseDatas(:,1:2); % epsilon in the model in section 5;
        %multiple = baseDatas(:,3); % determines how multiplicity is resolved;

        % X denotes the market type indicator;
        %X = zeros(size(baseDatas(:,4)));
        %for j=1:dimSX
        %    X = X + (j-1)*(baseDatas(:,4)>=(j-1)/dimSX).*(baseDatas(:,4)<j/dimSX) ;
        %end
        %betas0 = [0,theta0(3:end)]; % vector indicates beta_q for q=1,...,d_X

        % Initialize matrices that will contain both entry decision and market type
    %     dataP11 = NaN(size(baseDatas,1),size(PX,1));
    %     dataP10 = NaN(size(baseDatas,1),size(PX,1));
    %     for j=1:dimSX
    % 
    %         % Entry decision that indices {A_1=1,A_2=1}
    %         Entry11_aux = (epsilons(:,1) > theta0(1)-betas0(j)).*(epsilons(:,2) >  theta0(2)-betas0(j)) ;
    %         % should have mean = (1-(theta0(1)-betas0(j)))*(1-(theta0(2)-betas0(j)));
    % 
    %         % Entry decision that indices {A_1=1,A_2=0}
    %         Entry10_aux = (epsilons(:,1) > theta0(1)-betas0(j)).*(epsilons(:,2) <= theta0(2)-betas0(j)) + (epsilons(:,1) <= theta0(1)-betas0(j)).*(epsilons(:,2) <= theta0(2)-betas0(j)).*(multiple>prob10);
    %         % should have mean = (theta0(1)-betas0(j))*(theta0(2)-betas0(j))*(1-prob10) + (1-(theta0(1)-betas0(j)))*(theta0(2)-betas0(j));
    % 
    %         % entry (positive or zero) and market type (integer) information (It assumes P(X=x) is known)
    %         dataP11(:,j) = Entry11_aux.*(X == j-1)/PX(j); % indicates {A_1=1,A_2=1} and {X=j-1}
    %         dataP10(:,j) = Entry10_aux.*(X == j-1)/PX(j); % indicates {A_1=1,A_2=0} and {X=j-1}
    %     end

        % collect entry data in a single matrix;
        data =  W; %data = [dataP11,dataP10];
        dataP11 = W(:,1:dimSX);
        dataP10 = W(:,dimSX+1:end);

        % determine the value of theta of interest according to H0;
        theta_H0 = theta_H0_list(theta_index);

        % determines the value of the other thetas;
        if coordinate==1
            theta_other0 = theta0(2:end);
        else % coordinate = 2
            theta_other0 = [theta0(1),theta0(3:end)];
        end

        % indicates if H0 is true or not;
        %theta_in_IS = trueH0(theta_index);

        %% Compute test statistic for all tests, profile test statistic.

        % Choose starting values for minimization;
        starting_values = NaN(2,size(theta_other0,2)); % initialize;

        % first starting value: true parameter value (unfeasible in practice)
        starting_values(1,:) = theta_other0';

        % second starting value: good guess based on the theory (feasible in practice)
        starting_values(2,1) = 1 - mean(dataP11(:,1))/(1-theta_H0);
        % the remaining starting values are randomly chosen;
        for a=2:size(starting_values,2)
            if theta_other0(a)>0
                starting_values(2,a) = -(1/2)*(2-theta_H0-starting_values(2,1))+(1/2)*sqrt( (starting_values(2,1) - theta_H0)^2 + 4*mean(dataP11(:,a)) );
            else
                starting_values(2,a)=0;
            end
        end

        % Pick a large value as initial function value for the minimization
        min_value = 10^10;

        % for minimization we need the following constrains: Aineq*theta0?Bineq and (0,...,0)<=theta0<=(1,...1)
        Aineq = zeros(2*(size(theta0,2)-2),size(theta0,2)-1);
        Aineq(1:size(theta0,2)-2,2:size(theta0,2)-1) = eye(size(Aineq(1:size(theta0,2)-2,2:size(theta0,2)-1)));
        Aineq(size(theta0,2)-1:end,2:size(theta0,2)-1) = eye(size(Aineq(size(theta0,2)-1:end,2:size(theta0,2)-1)));
        Aineq(size(theta0,2)-1:end,1) = -1*ones(size(Aineq(size(theta0,2)-1:end,1)));
        Bineq = zeros(2*(size(theta0,2)-2),1);
        Bineq(1:(size(theta0,2)-2)) = theta_H0 * ones(size(Bineq(1:(size(theta0,2)-2))));

        min_outcomes = []; % This matrix will collect results of minimization;

        % solve numerical minimization for all starting values;
        for s=1:size(starting_values,1)
            [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_function(x,theta_H0,coordinate,data,p),starting_values(s,:),Aineq,Bineq,[],[],zeros(size(theta_other0)),ones(size(theta_other0)),[],options);
            % check whether minimization is successful and reduced value
            if Qn_aux  < min_value && bandera >= 1
                Qn_minimizer = theta_aux;
                min_value = Qn_aux;
            end
            % if minimization is successful, collect minimizer and its value;
            if bandera >= 1
                min_outcomes = [min_outcomes;[min_value,Qn_minimizer]]; %#ok<AGROW>
            end
        end
        % at the end of the minimizations, we should have a minimizer
        minQn = min_value;

        % Collect minimizers to estimate the indetified set (used for DR)
        min_outcomes = uniquetol2(min_outcomes,tolerance,'rows'); % set of minimizers;
        In_identified_set_hat = min_outcomes(min_outcomes(:,1)<=min(min_outcomes(:,1)) + tolerance , 2:end); % estimator of id set

        %% BP, i.e., projection.
        % By definition, this method requires gridsearch. It is feasible in simpler problems but not in this model and using Matlab.

        %% Subsampling test
    %     for bn_type=1:size(bn_list,1)
    % 
    %         for r=1:R
    %             % collect relevant data for SS minimization;
    %             data_SS = data(W3_SS(r,1:bn_list(bn_type)),:);
    % 
    %             % Starting values for SS minimization;
    %             starting_values =  uniquetol2([starting_values;Qn_minimizer],tolerance,'rows');
    % 
    %             % Pick a large value as initial function value for the minimization;
    %             min_value_SS = 10^10;
    % 
    %             % solve numerical minimization for all starting values;
    %             for s=1:size(starting_values,1)
    %                 [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_function(x,theta_H0,coordinate,data_SS,p),starting_values(s,:),Aineq,Bineq,[],[],zeros(size(theta_other0)),ones(size(theta_other0)),[],options);
    %                 % check whether minimization is successful and reduced value
    %                 if Qn_aux  < min_value_SS && bandera >= 1
    %                     minimizer = theta_aux;
    %                     min_value_SS = Qn_aux;
    %                 end
    %             end
    %             minQn_SS(bn_type,r) = min_value_SS;
    %         end
    % 
    %         % Compute SS critical value
    %         cn_SS(bn_type) = quantile(minQn_SS(bn_type,:),1-alpha);
    %     end
    %     % SS Test: Accept decision;
    %     accept_test_SS = repmat(minQn,1,size(bn_list,1)) <= cn_SS;

        % Starting values for PR minimization;
        starting_values =  uniquetol2([starting_values;Qn_minimizer],tolerance,'rows');

        %% MR test (and also DR and PR)
        for kappa_type =1:size(kappa_list,1)
            parfor r=1:R

                %- Step 1: Discard resampling method

                % Pick a large value as initial function value for the minimization;
                min_value_DR = 10^10;

                % Minimize but restricted to points of the estimated indentified set;
                for IDsetHat_index = 1:size(In_identified_set_hat,1)
                    minQn_DR_aux = min(min_value_DR,Qn_MR_function(In_identified_set_hat(IDsetHat_index,:),theta_H0,coordinate,data,kappa_list(kappa_type),p,k,W2_AA(r,:),1));
                end

                % compute simulated DR criterion function
                minQn_DR(kappa_type,r) = minQn_DR_aux;

                %- Step 2: Penalized resampling method;
                % Pick a large value as initial function value for the minimization;
                min_value_PR = 10^10;
                % check whether minimization is successful and reduced value
                for s=1:size(starting_values,1)
                    [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_MR_function(x,theta_H0,coordinate,data,kappa_list(kappa_type),p,k,W2_AA(r,:),2),starting_values(s,:),Aineq,Bineq,[],[],zeros(size(theta_other0)),ones(size(theta_other0)),[],options);
                    if Qn_aux  < min_value_PR && bandera >= 1
                        minimizer = theta_aux;
                        min_value_PR = Qn_aux;
                    end
                end

                % compute simulated PR criterion function
                minQn_PR(kappa_type,r) = min_value_PR ;

                %- Step 3: combine PR and DR to get MR
                minQn_MR(kappa_type,r) = min(minQn_DR(kappa_type,r),minQn_PR(kappa_type,r));
            end
            % Compute critical values
            cn_DR(kappa_type) = quantile(minQn_DR(kappa_type,:),1-alpha);
            cn_PR(kappa_type) = quantile(minQn_PR(kappa_type,:),1-alpha);
            cn_MR(kappa_type) = quantile(minQn_MR(kappa_type,:),1-alpha);
        end

        % DR, PR, and MR Tests: Accept decisions
        accept_test_DR = repmat(minQn,1,size(kappa_list,1)) <= cn_DR;
        accept_test_PR = repmat(minQn,1,size(kappa_list,1)) <= cn_PR;
        accept_test_MR = repmat(minQn,1,size(kappa_list,1)) <= cn_MR;

        %% Finish simulation: save and display time;
        %tElapsed = toc;

        % Save newly found results
        %cd results
        % load the corresponding matrices;
        %load results
        % replace the NaN with the corresponding results;
        %matrix_reject_test_SS(MC_index,theta_index,:) = ones(size(accept_test_SS)) - accept_test_SS; %#ok<SAGROW>
        matrix_reject_test_DR(1,theta_index) = ones(size(accept_test_DR)) - accept_test_DR; %#ok<SAGROW>
        matrix_reject_test_PR(1,theta_index) = ones(size(accept_test_PR)) - accept_test_PR; %#ok<SAGROW>
        matrix_reject_test_MR(1,theta_index) = ones(size(accept_test_MR)) - accept_test_MR;
        %save results matrix_reject_test*
        %cd ..

        % Display overall time of the MCrun in minutes
        %if display_info == 1
        %    disp(['- MCrun done in ' num2str(tElapsed/60) ' mins.']);
        %end
        %fprintf('iter %9.4f of %9.4f done \n', [theta_index,size(theta_H0_list,1)])

        % If accept MR test, then we have found the lower bound.
        if matrix_reject_test_MR(1,theta_index) == 0
            break;
        end
    end
 
    ind = find(isnan(matrix_reject_test_MR));
    matrix_reject_test_MR(1,ind) = 0;
    
    BCS_output.MRreject_ind = matrix_reject_test_MR;
    BCS_output.grid = theta_H0_list;
    BCS_output.MRaccept = theta_H0_list(find(matrix_reject_test_MR ==0)).'; 
    BCS_CI =[min(BCS_output.MRaccept),max(BCS_output.MRaccept)];
    BCS_output.flag_conv = 1;
catch
    % Sometimes get an error in fmincon due to numerical problems.  Output
    % NaN for this MC and flag.
    BCS_CI = [nan,nan];
    BCS_output.flag_conv = 0;
end




end