%This Matlab function implements the subvector CC test of Cox and Shi
%(2020). H_0: there exists delta such that C*delta<=m

%Inputs: 
%C: estimator of C
%m: estimator of m
%V: a consistent estimator for the (asymptotic) variance of m. This can be
%   a conditional variance estimator.
%alpha: significance level

%Outputs:
%T_CC: the quasilikelihood ratio statistic
%cv_CC: the conditional chi-squared critical value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T_CC,cv_CC, dof_n] = func_subCC(C, m, V,alpha)
  
        d_nuis = size(C,2); %number of nuisance parameters
        d_ineq = length(m); %number of inequalities
        invVar = inv(V);
        
        %Perform CC Test
        
        %Finding T (Test Statistic)
        H = 2*[zeros(d_nuis,d_ineq+d_nuis);zeros(d_ineq,d_nuis),invVar];
        f = -2*[zeros(d_nuis,1);invVar*m];
        
        Aaug = [C,-eye(d_ineq)];
        
        options = optimset('TolCon',1e-18,'MaxIter',500,'display','off');
        lp_options = optimset('display','off','algorithm','dual-simplex');
%        lp2_options = optimset('algorithm','interior-point','maxiter',2000,'TolCon',1e-4,'display','off');
        
        [y_star,~] = quadprog(H,f,Aaug,zeros(d_ineq,1),[],[],[],[],[],options);
 %       delta_star = y_star(1:d_nuis);
        m_star = y_star(d_nuis+1:d_nuis+d_ineq);
        
        T_CC = (m-m_star)'*invVar*(m-m_star);  %CC Test Statistic
        
        %Find the Degree of Freedom:
        
        %The dof is equal to the maximum number of linearly independent vectors
        %in {b: b>=0,C'*b=0, m_star'*b=0,1'b=1}.
        %We use the enumerating method to find all the implicit equalities
        A_ineq0 = - eye(d_ineq);
        b_ineq0 = zeros(d_ineq,1);
        A_eq0 = [C';m_star';ones(1,d_ineq)];
        A_eq = [C';ones(1,d_ineq)];
%        b_eq0 = [zeros(d_nuis,1);0;1];
        b_eq = [zeros(d_nuis,1);1];
        
        %By definition, m_star'b>=0 as long as b>=0,C'b=0 and 1'b=1. If we
        %find that m_star'b>0 for all such b's, there must be no active inequalities among
        %B(C)'m>=0. Thus, we use the folliwng LP to identify this case.
        
    
        [~,Vmu_min,flag] = linprog(m_star,A_ineq0,b_ineq0,A_eq,b_eq,[],[],lp_options);

        if (Vmu_min>=0.00005)
            dof_n=0;
        elseif (flag==-2)
            dof_n = 0;
        else
            %If it is found that dof_n>0, proceed to find all the implicit
            %equalities:
            
            nb_ineq_min = NaN(length(b_ineq0),1);
            counter=1;
            for bj = 1:d_ineq   %Go through each j one by one
           %find the largest b_j allowed by b>=0,C'b=0, 1'b=1,m_star'b=0:

                [~,nbj_min] = linprog(A_ineq0(bj,:)',[],[],...
                    A_eq0,[zeros(d_nuis,1);Vmu_min;1],...
                    b_ineq0-(1e-7),b_ineq0+1,lp_options);
                %nbj_min stands for minumum of negative bj allowed (if it
                %is close to zero, the jth ineq is an implicit equality).

                if isempty(nbj_min)
                    nb_ineq_min(bj) = b_ineq0(bj)-1;
                    counter=counter+1;
                    %if no feasible solution, consider bj a strict inequality to be on the safe side.
                else
                    nb_ineq_min(bj) = nbj_min;
                end
            end
                       
            
            A_impeq = A_ineq0((b_ineq0 - nb_ineq_min)<(1e-4)/2,:);
            %Collect the rows of A_ineq corresponding to all the implicit
            %equalities.
            
            A_full_eq = [A_impeq;A_eq0];
            rkA = rank(A_full_eq,(1e-4)/2);  %The tolerance in the rank 
                                             %function is important sometimes 
                                             %because the calculation above may
                                             %not be previse enough
            
            dof_n = d_ineq - rkA + 1;
        end
        cv_CC = chi2inv(1-alpha,dof_n); 

    end
    
 