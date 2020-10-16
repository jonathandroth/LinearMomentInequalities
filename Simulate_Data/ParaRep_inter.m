%This finds the Parametric Representation of Polyhedrons

function [RA_ineq,Rb_ineq,Trans_A,Trans_b] = ParaRep_inter(A_eq0,b_eq0,A_ineq0,b_ineq0)

%Apply Gauss-Jordan Reduction if there are equalities
if min(size(A_eq0))>0

    [RA_ineq, Rb_ineq,Trans_A,Trans_b] = GaJo_inter(A_eq0,b_eq0, A_ineq0,b_ineq0);
else
    RA_ineq = A_ineq0;
    Rb_ineq = b_ineq0;
    Trans_A = eye(size(A_ineq0,2));
    Trans_b = zeros(size(A_ineq0,2));
end

no_eq_ind = 0;
while no_eq_ind == 0 && min(size(RA_ineq))>0

             
        %Step 3: Check for implicit equalities by solving
       
        rAA = size(RA_ineq,1);
        cAA = size(RA_ineq,2);
        options = optimset('display','off');
        [lam,f,exitflag] = linprog(Rb_ineq,[],...
            [],[RA_ineq';ones(1,rAA)],[zeros(cAA,1);1],zeros(rAA,1),[],options);
        
        if exitflag == 1 && abs(f)<1e-10
            A_ineq0 = RA_ineq((lam==0),:);
            b_ineq0 = Rb_ineq((lam==0));
            A_eq0 = RA_ineq((lam>0),:);
            b_eq0 = Rb_ineq((lam>0));

            [RA_ineq, Rb_ineq,Trans0_A,Trans0_b] = ...
                GaJo_inter(A_eq0, b_eq0,A_ineq0,b_ineq0);
            Trans_A = Trans_A*Trans0_A; 
            Trans_b = Trans_A*Trans0_b + Trans_b;
        else
            no_eq_ind=1;
        end
end
    
end