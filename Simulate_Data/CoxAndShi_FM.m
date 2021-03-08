function [CS_test, QLR_stat, CS_crit] = CoxAndShi_FM(Y,X,Sigma,alpha)

%Find the Fourier Motzkin matrices implied by the systes E[Y-X\delta] <=0
[A_fm, b_fm] = findFMMatAndVec(X);

%Normalize A_fm and b_fm so each row of A has norm of 1 
    %This helps with numerica precision
row_norms = arrayfun(@(n) norm(A_fm(n,:)), 1:size(A_fm,1));
A_fm = A_fm ./ repmat(row_norms',1, size(A_fm,2));
b_fm = b_fm ./ row_norms';

%Run the unconditional Cox and Shi test using the implied matrices
[CS_test, QLR_stat, CS_crit] = CoxAndShi_NoNuisance(Y,A_fm,b_fm,Sigma,alpha);

end
