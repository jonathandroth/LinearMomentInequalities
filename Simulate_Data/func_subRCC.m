%This Matlab function implements the subvector RCC test of Cox and Shi
%(2020). H_0: there exists delta such that C*delta<=m

%Inputs: 
%C:    estimator of C
%m:    estimator of m
%Var:  a consistent estimator for the variance of m (not for sqrt(n)*m). 
%       This can be a conditional variance estimator.
%alpha: significance level

%Outputs:
%T_RCC: the quasilikelihood ratio statistic
%cv_RCC: the refined conditional chi-squared critical value
%        only apply the refinement when dof_n=1 and T_RCC \in
%        [chi2inv(1-2*alpha,1),chi2inv(1-alpha,1)]
%cv_CC: the conditional chi-squared critical value (unrefined)
%cv_CC2: the conditional chi-squared critical value computed after vertices
%        enumeration.
%Edim:   dimension of the vertice enumeration result.

%This code calls for the user-defined functions ParaRep_inter() and
%inequalities2vertices()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T_RCC,cv_RCC,cv_CC,cv_CC2,Edim, dof_n] = func_subRCC(C, m, Var,alpha)
  
        d_nuis = size(C,2); %number of nuisance parameters
        d_ineq = length(m); %number of inequalities
        invVar = inv(Var);
        
        %Perform CC Test first
        
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
        
        T_RCC = (m-m_star)'*invVar*(m-m_star);  %CC Test Statistic also is the RCC test statistic
        
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

        
        %% The Refinement
        
        if dof_n==1 && T_RCC<=chi2inv(1-alpha,1) && T_RCC> chi2inv(1-2*alpha,1)
            %This is the case where the refinement could make a difference
            
            %We enumerate the vertices of the polyhedron: {eta:eta>=0, C'*eta=0,1'eta=1}
            %The vertices form the rows of the E matrix such that E*m>=0 is the
            %implied set of moment inequalities.
            A_ineq = - eye(size(C,1));
            A_eq = [C';ones(1,size(C,1))];
            b_eq = [zeros(size(C,2),1);1];
            b_ineq = zeros(size(C,1),1);
            
            %Finding the parametric form of the polyhedron so that
            %inequalities2vertices() can be applied on the full-dimensional
            %projection:
            [RA_ineq,Rb_ineq,Trans_mat,Trans_b] = ParaRep_inter(A_eq, b_eq,A_ineq,b_ineq);
            [V,~] = inequalities2vertices(RA_ineq,Rb_ineq);
            Vfull = V*Trans_mat' + repmat(Trans_b',size(V,1),1); %transform back from the full-dimensional projection
            Edim = [size(Vfull,1),size(Vfull,2)];
            
            %CC Critical Value
            Vmu = Vfull*m;
            dof_n2 = rank(Vfull(Vmu<(1e-4)/2,:));
            cv_CC2 = chi2inv(1-alpha,dof_n2);
            
            %RCC critical Value
                [~,jmin] = min(Vmu);
                a_jmin = Vfull(jmin,:)';
                zjX = NaN(length(Vmu),1);
                for j=1:length(Vmu)
                    a_j = Vfull(j,:)';
                    if abs(sqrt(a_jmin'*Var*a_jmin*(a_j'*Var*a_j))...
                                         -a_jmin'*Var*a_j) < (1e-4)/2
                        zjX(j)=inf;
                    else
                        zjX(j) = sqrt(a_jmin'*Var*a_jmin)*(a_j'*m)/...
                            (sqrt(a_jmin'*Var*a_jmin*(a_j'*Var*a_j))-a_jmin'*Var*a_j);
                    end
                end
                zX = max(0,min(zjX));
                beta = 2*alpha*normcdf(zX);
                cv_RCC = chi2inv(1-beta,1);
        else
            cv_RCC = cv_CC;
            Edim = [NaN,NaN];
            cv_CC2 = cv_CC;
        end

        
        

    end
