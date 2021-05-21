for ds = 1:numdatasets

ds
   X_T = X_T_cell{ds,1};
   y_T = y_T_cell{ds,1};
   Sigma = Sigma_conditional_cell{ds,1}; 
    

   %For the KMS/AMS and cond'l/hybrid, 
   %we do a change of basis such that M*theta yields l*theta in the first
   %row (I construct an orthonormal basis for the rest of the basis
   %vectors, although this is not strictly necessary)
   M = [l' ; null( repmat(l',size(l)) )'];
   X_T = X_T * M^(-1);
   
   
   %Do AS and KMS for specs with <9 parameters
   if size(X_T,2) < 9
   try
       [as_ci,as_output] = projected_AS_or_KMS(y_T, X_T, 500, Sigma,[1;zeros(size(X_T,2)-1,1)], NaN, 'AS');
       
       %If have a convergence error for either bound, set to NaN
       if (as_output.flagL_EAM ~= 1) || (as_output.flagU_EAM ~= 1)
           as_ci = [NaN, NaN];
       end
       
   catch
   end
   
   try
   
       [kms_ci,kms_output] = projected_AS_or_KMS(y_T, X_T, 500, Sigma,[1;zeros(size(X_T,2)-1,1)], NaN, 'KMS');
       
       %If have a convergence error for either bound, set to NaN
       if (kms_output.flagL_EAM ~= 1) || (kms_output.flagU_EAM ~= 1)
           kms_ci = [NaN, NaN];
       end
       
   catch
   end
   
   if(kms_ci(1) < as_ci(1) || kms_ci(2) > as_ci(2) )
       warning('AS bound tighter than KMS bound')
   
   end
   
   end
end