function [y_T,X_T,Sigma] = add_mean_to_moments_fn(c, y_T,X_T,Sigma)


%c = 0.00001;
   
   %Add c * the mean to each of the moments
   k = size(X_T,1);
   A = ( eye(k,k) + c*ones(k,k) );
   
   X_T = A * X_T;
   y_T = A * y_T;
   Sigma = A * Sigma * A';
   
   %Renormalize the moments to have variance one
   D_minushalf = diag( diag(Sigma).^(-1/2) );
   X_T = D_minushalf * X_T;
   y_T = D_minushalf * y_T;
   Sigma = D_minushalf * Sigma * D_minushalf;