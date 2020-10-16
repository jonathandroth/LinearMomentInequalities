%This returns the Gauss-Jordan Reduction of a Set of Linear Inequalities
%and Equalities

%This should only be used if rank(A_eq) = rank([A_eq,b_eq])

function [RA_ineq,Rb_ineq,Trans_mat,Trans_b] = GaJo_inter(A_eq, b_eq,A_ineq,b_ineq)

if rank([A_eq,b_eq])-rank(A_eq)==0
[RAb_eq,pAb_eq] = rref([A_eq,b_eq]);   %Reduced row Echelon Form of [A_eq,b_eq]

rL1 = sum(abs(RAb_eq),2);
RAb_eq(rL1==0,:)=[];  %remove rows that are all zeros

Rb_eq = RAb_eq(:,size(RAb_eq,2));  %transformed b
                                   
RA_eq_npiv = RAb_eq(:,1:size(RAb_eq,2)-1);
RA_eq_npiv(:,pAb_eq) = []; %Nonpivotal Columns of RAb_eq
                          %x_piv = -RA_eq_npiv*x_npiv + Rb_eq;                         
A_ineq_piv = A_ineq(:,pAb_eq); %pivotal columns of RA_ineq
A_ineq_npiv = A_ineq;
A_ineq_npiv(:,pAb_eq) = [];    %Nonpivotal columns of A_ineq
                       %A_ineq_npiv*x_npiv + A_ineq_piv*x_piv <= b_ineq.
%Thus, (A_ineq_npiv - A_ineq_piv*RA_eq_npiv)*x_npiv<=b_ineq - A_ineq_piv*Rb_eq

RA_ineq = A_ineq_npiv - A_ineq_piv*RA_eq_npiv;

Rb_ineq = b_ineq - A_ineq_piv*Rb_eq;


%Matrix to transform the reduced variables to original vector of variables.
ind = 1:size(A_eq,2);
npA_eq = ind;
npA_eq(pAb_eq)=[];   %indices for non-pivotal columns

Trans_mat = NaN(size(A_eq,2),length(npA_eq));

Trans_mat(pAb_eq,:) = -RA_eq_npiv;  %matrix to write full vector of variables
                          %in terms of not eliminated ones
Trans_mat(npA_eq,:) = eye(length(npA_eq));    

Trans_b = NaN(size(A_eq,2),1);
Trans_b(pAb_eq) = Rb_eq;
Trans_b(npA_eq) = zeros*npA_eq;

%x = Trans_mat*x_npiv  + Trans_b 
else
    error('Inequality System Has No Solution')
end

                       
                       
                       
                       
                       
                       
                       
                       
                       
                       
                       