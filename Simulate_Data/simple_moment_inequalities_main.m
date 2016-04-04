

load( '../../Output/Simulated_Data/ds1.mat' )

T = size( Pi_array ,3);
J = size( Pi_array, 2);
F = size( Pi_array, 1);


T_array = repmat( (1:T), F*J ,1  );
T_array = T_array(:);
T_array = reshape( T_array, F, J,T);

%Create an array of profits for the weight class one weight ahead (behind). Set to 0
%if the maximum (minimum) weight class
Pi_gplus1_array = [ Pi_array(:,2:J,:), zeros(F,1,T)];
Pi_gminus1_array = [ zeros(F,1,T), Pi_array(:,1:(J-1),:)];


%Put into "long" format
Pi_vec = Pi_array(:);
Pi_gplus1_vec = Pi_gplus1_array(:);
Pi_gminus1_vec = Pi_gminus1_array(:);

J_t_vec = J_t_array(:);
J_tminus1_vec =  J_tminus1_array(:);
G_vec = G_array(:);
F_vec = F_array(:);
T_vec = T_array(:);

%Specify number of moments
M= 6;

N = size(J_t_vec,1); % compute number of obs

%Create C_jft mat
C_mat = NaN( N , M);

C_mat(:,1) = (J_t_vec == 1 & J_tminus1_vec == 1);
C_mat(:,2) = (J_t_vec == 0 & J_tminus1_vec == 1);
C_mat(:,3) = (J_t_vec == 1 & J_tminus1_vec == 0);
C_mat(:,4) = (J_t_vec == 0 & J_tminus1_vec == 0);
C_mat(:,5) = (J_t_vec == 1 & J_tminus1_vec == 0);
C_mat(:,6) = (J_t_vec == 1 & J_tminus1_vec == 0);

moment_has_lambda = [1,1,0,0,0,0];%vector indicating whether the moment "needs a lambda"

% Create a matrix where the i,j th entry is the value of delta-pi for observation
% i if C_ij == 1, and is 0 otherwise

%The deltapi is just Pi_vec for the first 4 moments. For the 5th moment, this is
%(Pi_vec - Pi_gminus1_vec); for the 6th, this is (Pi_vec - Pi_gplus1_vec)

Pi_cond_mat = [repmat(Pi_vec, 1, M-2), Pi_vec - Pi_gminus1_vec, Pi_vec - Pi_gplus1_vec] .* C_mat;

% Create a matrix where the i,j th entry is the value that multiplies theta_g in the moment for observation
% i if C_ij == 1, and is 0 otherwise

%For moments 1:4, the value that multiplies theta_g is +g if the product is
%in, and -g if it is out. For the last two, it is +1 and -1
G_cond_mat = repmat(G_vec, 1, M) .* C_mat;


%Avg the values of Pi and G for each (t,g) group. All of the moments will
%be of the form: M = Pi_for_moments - constant1 - constant2 * G_for_moments
Grouping_mat = [ T_vec, G_vec];

% tic
% for m = 1:M
%     Pi_for_moments = accumarray( grp2idx( num2str(Grouping_mat)), Pi_cond_mat(:,m), [], @mean) ;
%     G_for_moments = accumarray( grp2idx( num2str(Grouping_mat)), G_cond_mat(:,m), [], @mean) ;
% end
% toc

Pi_for_moments = grpstats( Pi_cond_mat, {T_vec, G_vec}, @mean);
G_for_moments = grpstats( G_cond_mat, {T_vec, G_vec}, @mean);
Const_for_moments = grpstats( C_mat, {T_vec, G_vec}, @mean);

%Create a matrix that indicates whether each moment should be multiplied by
%lambda or not
%This is a matrix of size #grps x M where each column is all 1s if the mth moment
%relates to a case where J_tminus1= 1 and is 0 otherwise

Lambda_indicator_mat = repmat( moment_has_lambda , size(Pi_for_moments,1) , 1);

%Create a matrix that has a column of 1s for moments in which (pi - fc) >0 
% and has a columns of -1s for moments in which (pi -fc) <0
Sign_moment_mat = repmat( [1,-1,1,-1,1,1], size(Pi_for_moments,1) ,1 );

%Create a  fn that returns a (F*J) x M matrix where the (i,j)th entry is the jth moment for
%grp i

%NEED TO FIX LAST TWO MOMENT INEQUALITIES. SHOULD HAVE ONES INSTEAD OF
%THETA_G

%INSTEAD OF DOING SIGN_MOMENTS_MAT, CAN HAVE CONST_FOR_MOMENTS AND
%G_FOR_MOMENTS HAVE A DIFFERENT SIGN

%Also need to fix the grpstats so that they avg only over the things that
%are positive (i.e. the observations for which C > 0)
Moments_mat_fn = @(theta_c, theta_g, lambda) (Pi_for_moments - ( lambda * Lambda_indicator_mat ...
    + (1 - Lambda_indicator_mat) ) .* ( theta_c * Const_for_moments + theta_g * G_for_moments)) ...
    .* Sign_moment_mat;



Moments_mat_true = Moments_mat_fn(129.73, -21.38, 0.386);

sum( Moments_mat_true)

%  testJ = (Pi_array -  (J_tminus1_array * lambda + (1 - J_tminus1_array))...
%      .* (theta_c + G_array *theta_g)      ) >= 0 ;
%  
%  diffJ = testJ - J_t_array;
%  sum( abs( diffJ(:) ) )
% 
