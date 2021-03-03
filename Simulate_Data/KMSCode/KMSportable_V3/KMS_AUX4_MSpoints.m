function [Mpoints] = KMS_AUX4_MSpoints(theta_Astep)
% Receive points used in A-step
X = theta_Astep;
[Neval,dim_p] = size(X);

% Create a tentative matrix of Neval x dim_p x dim_p+1
ms_points = zeros(Neval,dim_p,dim_p+1);

% Create initial start points for EI maximization
for k=1:dim_p    
[temp,I] = sortrows(X,k);   % sort k-th column 
ms_points(:,:,k) = temp;    % for now fill in the sorted eval points

% Replace k-th coordinate with averages with nearest neighbors.
% See Shonlau, Welch, and Jones (98) page 16. (Points A and B)
temp1 = [(temp(1,k)+temp(2,k))/2; (temp(2:end,k)+ temp(1:end-1,k))/2];
temp2 = [temp1(2:end);(temp(end-1,k)+temp(end,k))/2];
sel = (rand(Neval,1)>0.5);

% Randomize between two candidates
ms_points(:,k,k) = sel.*temp1 + (1-sel).*temp2; 

% Use the ones that were not selected for dim_p+1-th evaluation point (Point E)
temp3 = (1-sel).*temp1 + sel.*temp2;
ms_points(I,k,dim_p+1) = temp3;
end

% Reshape the evaluation points as (Neval*dim_p+1) x dim_p
for k=1:dim_p+1
Mpoints((k-1)*Neval+1:k*Neval,:)=ms_points(:,:,k);
end
end

