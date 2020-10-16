function [T_CC, cv_RCC, cv_CC] = rcc_test_fn(m_n, B_z, C_z, d_z, Sigma, n, CConly, alpha)

%XX need to either define T_CC_mat here or get rid of the sim
%subsetting

%Extract dimensions from the arguments;
k = length(d_z);
d_m = size(B_z,2);
p = size(C_z,2);

invCondVar = inv(Sigma);

%The moments have the term C_z delta, so since C_z is k x p, delta is px1
%Sigma is the variance of m and hence is d_m x d_m
%It appears that d_x+1 corresponds with p and 2*dzf with d_m
%It appears that g_nbar is the m_n value, and d_Z =0 in their code
% It appears that X_nbar corresponds with -C_z? 

%Finding T (Test Statistic)
%H = 2*[zeros(d_x+1,d_x+1+2*dzf);zeros(2*dzf,d_x+1),invCondVar];
%f = -2*[zeros(p,1);invCondVar*g_nbar];
%Aaug = [X_nbar,-eye(2*dzf)];

%The optimization variable is the stacked vector (delta',mu')'
%Quadprog(H,f,A,b) computes min_x 1/2 * x'Hx + f'x s.t. A x <= b

%We construct H to be the block matrix with 0s corresponding with the delta
%component of x and 2*invCondVar in the bottom right block corresponding with mu
H = 2*[zeros(p, p + d_m); zeros(d_m, p), invCondVar];
%We construct f to have 0s corresponding with delta and -2* invCondvar *
%m_n corresponding with mu
f = -2*[zeros(p,1);invCondVar*m_n];
%We construct A to specify the constraint (-C_z delta + B_z * mu) <= d_Z  
Aaug = [-C_z,B_z];
b = d_z;
options = optimset('TolCon',1e-18,'MaxIter',500,'display','off');
lp_options = optimset('display','off','algorithm','dual-simplex');

% [y_star,~] = quadprog(H,f,Aaug,zeros(dzf*2,1),[],[],[],[],[],options);
% delta_star = y_star(1:d_x+1);
% mu_star = y_star(d_x+2:d_x+1+2*dzf);
% T_CC_mat(sim) = n*(g_nbar-mu_star)'*invCondVar*(g_nbar-mu_star);  %CC Test Statistic

[y_star,~] = quadprog(H,f,Aaug,b,[],[],[],[],[],options);
delta_star = y_star(1:p);
mu_star = y_star(p+1:p+d_m);
T_CC = n*(m_n-mu_star)'*invCondVar*(m_n-mu_star);  %CC Test Statistic


        %Find the Degree of Freedom:
        
        %The dof is equal to the maximum number of linearly independent vectors
        %in {b: b>=0,X_nbar'*b=0, mu_star'*b=0,1'b=1}.
        %We use the enumerating method to find all the implicit equalities
        A_ineq0 = - eye(size(-C_z,1));
        b_ineq0 = zeros(size(-C_z,1),1);
        A_eq0 = [-C_z';mu_star';ones(1,size(-C_z,1))];
        A_eq = [-C_z';ones(1,size(-C_z,1))];
        b_eq0 = [zeros(size(-C_z,2),1);0;1];
        b_eq = [zeros(size(-C_z,2),1);1];
        
        [~,Vmu_min,exit] = linprog(mu_star,A_ineq0,b_ineq0,A_eq,b_eq,[],[],lp_options);
        if Vmu_min>=0.00005
            dof_n = 0;
        else
            b_ineq_min = NaN(length(b_ineq0),1);
            counter=1;
            for bj = 1:length(b_ineq0)
                [~,bj_min] = linprog(A_ineq0(bj,:)',[],[],...
                    A_eq0,[zeros(size(-C_z,2),1);Vmu_min;1],b_ineq0-(1e-7),b_ineq0+1,lp_options);
                
                if isempty(bj_min) % if no feasible solution, consider bj a strict inequality to be on the safe side.
                    b_ineq_min(bj) = b_ineq0(bj)-1;
                    counter=counter+1;
                else
                    b_ineq_min(bj) = bj_min;
                end
            end
%            counter_mat(sim) = counter;
            
            
            A_impeq = A_ineq0((b_ineq0 - b_ineq_min)<(1e-4)/2,:);
            
            A_full_eq = [A_impeq;A_eq0];
            rkA = rank(A_full_eq);
            
            dof_n = size(A_full_eq,2)-rkA+1;
        end
        cv_CC = chi2inv(1-alpha,dof_n);
        
        %Compute the B matrix and refinement when refinement might make a difference.
        if CConly==0 && dof_n==1 && T_CC<=chi2inv(1-alpha,1) &&...
                T_CC> chi2inv(1-2*alpha,1)
            %We enumerate the vertices of the polyhedron: {b:b>=0, X_nbar'*b=0,1'b=1}
            %The vertices form the B matrix such that B*mu_star>=0 is the
            %implied set of moment inequalities.
            A_ineq = - eye(size(-C_z,1));
            A_eq = [-C_z';ones(1,size(-C_z,1))];
            b_eq = [zeros(size(-C_z,2),1);1];
            b_ineq = zeros(size(-C_z,1),1);
            
            %Finding the parametric form of the polyhedron so that
            %inequalities2vertices() can be applied on the full-dimensional
            %projection:
            [RA_ineq,Rb_ineq,Trans_mat,Trans_b] = ParaRep_inter(A_eq, b_eq,A_ineq,b_ineq);
            [V,nr] = inequalities2vertices(RA_ineq,Rb_ineq);
            Vfull = V*Trans_mat' + repmat(Trans_b',size(V,1),1); %transform back from the full-dimensional projection
            Bdim = [size(Vfull,1),size(Vfull,2)];
            
            %CC Critical Value
            Vmu = Vfull*mu_star;
            dof_n = rank(Vfull(Vmu<(1e-4)/2,:));
            cv_CC2 = chi2inv(1-alpha,dof_n);
            
            %RCC critical Value
                [~,jmin] = min(Vmu);
                a_jmin = Vfull(jmin,:)';
                zjX = NaN(length(Vmu),1);
                for j=1:length(Vmu)
                    a_j = Vfull(j,:)';
                    if abs(sqrt(a_jmin'*Sigma*a_jmin*(a_j'*Sigma*a_j))-a_jmin'*Sigma*a_j) < (1e-4)/2
                        zjX(j)=inf;
                    else
                        zjX(j) =sqrt(n)*sqrt(a_jmin'*Sigma*a_jmin)*(a_j'*mu_star)/...
                            (sqrt(a_jmin'*Sigma*a_jmin*(a_j'*Sigma*a_j))-a_jmin'*Sigma*a_j);
                    end
                end
                zX = max(0,min(zjX));
                beta = 2*alpha*normcdf(zX);
                cv_RCC = chi2inv(1-beta,1);
        else
            cv_RCC = cv_CC;
        end


end