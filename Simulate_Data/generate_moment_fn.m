
%This function takes the number of a dataset

%It returns as its output two functions of (theta_c, theta_g, lambda) that give
%the matrix of moments for each jxt cell. The first, Moments_mat_fn returns
%the matrix for the first 6 moments; the second Interacted_Moments_mat_fn,
%gives the first 6 momements as well as those interacted with eta+ and eta-

%Note: this function assumes the working directory is two levels below the
%main directory (e.g. in Code/something)

function [Moments_mat_fn, Interacted_Moments_mat_fn, Y, A_g, A_c, Y_basic, A_g_basic, A_c_basic]...
    = generate_moment_fn(F_array, G_array, Eta_shocks_array, Pi_array, J_t_array, J_tminus1_array)


T = size( Pi_array ,3);
J = size( Pi_array, 2);
F = size( Pi_array, 1);


T_array = repmat( (1:T), F*J ,1  );
T_array = T_array(:);
T_array = reshape( T_array, F, J,T);



%Construct the eta+ and eta- instrument
Eta_shocks_vec = Eta_shocks_array(:);
Eta_plus_vec = max(Eta_shocks_vec,0);
Eta_minus_vec = max(-Eta_shocks_vec,0);

%Put into "long" format
Pi_vec = Pi_array(:);

% %Create a J_array for g+1
% J_gp1_t_array = [ J_t_array(:,(2:size(J_t_array,2)),:) , zeros(size(J_t_array,1) , 1, size(J_t_array,3) ) ];
% J_gm1_t_array = [ zeros(size(J_t_array,1) , 1, size(J_t_array,3) ) , J_t_array(:,(1:(size(J_t_array,2)-1) ),:) ];
% J_gp1_tm1_array = [ J_tminus1_array(:,(2:size(J_tminus1_array,2)),:) , zeros(size(J_tminus1_array,1) , 1, size(J_tminus1_array,3) ) ];
% J_gm1_tm1_array = [ zeros(size(J_tminus1_array,1) , 1, size(J_tminus1_array,3) ) , J_tminus1_array(:,(1:(size(J_tminus1_array,2)-1) ),:) ];


J_t_vec = J_t_array(:);
J_tminus1_vec =  J_tminus1_array(:);

G_vec = G_array(:);

F_vec = F_array(:);
T_vec = T_array(:);





%%

%Specify number of moments
M= 6;

N = size(J_t_vec,1); % compute number of obs

%Create C_jft mat
C_mat = NaN( N , M);

C_mat(:,1) = (J_t_vec == 1 & J_tminus1_vec == 1);
C_mat(:,2) = (J_t_vec == 0 & J_tminus1_vec == 1);
C_mat(:,3) = (J_t_vec == 1 & J_tminus1_vec == 0);
C_mat(:,4) = (J_t_vec == 0 & J_tminus1_vec == 0);

moment_has_lambda = [1,1,0,0,0,0];%vector indicating whether the moment "needs a lambda"



%% Calculate profits for the plus-weight and minus-weight moments

%The allproducts f function takes an an argument a FxJxT panel
%It creates a matrix where each row is the row in the panel in   
allproducts_f_fn = @(panel) panel( sub2ind( size(panel), repmat(F_vec, 1, J), repmat( (1:J), size(F_vec,1), 1) ,repmat(T_vec, 1, J)) );

G_allf = allproducts_f_fn(G_array);
J_t_allf = allproducts_f_fn(J_t_array);
J_tminus1_allf = allproducts_f_fn(J_tminus1_array);
Pi_allf = allproducts_f_fn(Pi_array);

%Create indicators for the products that have weight higher and lower
upper_weights = G_allf > repmat(G_vec, 1, J) ;
lower_weights = G_allf < repmat(G_vec, 1, J) ;

%Create indicators for the products we want to sum over for the upper and
%lower weights
%These are products 
add_upper = upper_weights .* (J_t_allf == 0) .* (J_tminus1_allf == 0) ;
add_lower = lower_weights .* (J_t_allf == 0) .* (J_tminus1_allf == 0) ;

%Average pi over the eligible upper products (and set to 0 if no eligible
%upper products
Pi_upper = sum( Pi_allf .* add_upper  ,2) ./ max( 1, sum( add_upper, 2));
%Average pi over the elibible lowe products
Pi_lower = sum( Pi_allf .* add_lower  ,2) ./ max( 1, sum( add_lower, 2));

%Average G over the eligible upper products
G_upper = sum( G_allf .* add_upper  ,2) ./ max( 1, sum( add_upper, 2));
%Average G over the elibible lowe products
G_lower = sum( G_allf .* add_lower  ,2) ./ max( 1, sum( add_lower, 2));


%We are in condition 5 if we added product j and there is at least one
%eligible lower weight product
C_mat(:,5) = (J_t_vec == 1 & J_tminus1_vec == 0) & ( sum( add_lower, 2) > 0 );
C_mat(:,6) = (J_t_vec == 1 & J_tminus1_vec == 0) & ( sum( add_upper, 2) > 0 );

%%



% Create a matrix where the i,j th entry is the value of delta-pi for observation
% i if C_ij == 1, and is 0 otherwise

%For the first four moments, deltapi is pi if J_t =1 ,and -pi if J_t = 0
%the first 4 moments. For the 5th moment, this is Pi_vec - the average for the eligible lower products
%; for the 6th, this is pi_vec - the average_for the upper products

Pi_cond_mat = [Pi_vec, -Pi_vec, Pi_vec, -Pi_vec,...
    Pi_vec - Pi_lower, Pi_vec - Pi_upper] .* C_mat;

%Interact with eta_plus and eta_minus

Pi_cond_mat = [Pi_cond_mat,...
               Pi_cond_mat .* repmat( Eta_plus_vec, 1, size(Pi_cond_mat,2) ),...
               Pi_cond_mat .* repmat( Eta_minus_vec, 1, size(Pi_cond_mat,2) )];
% Create a matrix where the i,j th entry is the value that multiplies theta_g in the moment for observation
% i if C_ij == 1, and is 0 otherwise

%For moments 1:4, the value that multiplies theta_g is +g if the product is
%in, and -g if it is out. For the last two, it is -1 and +1
G_cond_mat = [-G_vec, G_vec, -G_vec, G_vec, G_lower - G_vec , G_upper - G_vec]  .* C_mat;

G_cond_mat = [G_cond_mat,...
               G_cond_mat .* repmat( Eta_plus_vec,1, size(G_cond_mat,2) ),...
               G_cond_mat .* repmat( Eta_minus_vec, 1, size(G_cond_mat,2) )];

% Create a matrix where the i,j th entry is the value that multiplies theta_c in the moment for observation
% i if C_ij == 1, and is 0 otherwise
           
Const_cond_mat = repmat([-1,1,-1,1,0,0], size(C_mat,1),1) .* C_mat;

Const_cond_mat = [Const_cond_mat,...
               Const_cond_mat .* repmat( Eta_plus_vec,1, size(Const_cond_mat,2) ),...
               Const_cond_mat .* repmat( Eta_minus_vec, 1, size(Const_cond_mat,2) )];



% %Avg the values of Pi and G for each (t,g) group. All of the moments will
% %be of the form: M = Pi_for_moments - constant1 - constant2 * G_for_moments
% Grouping_mat = [ T_vec, G_vec];

% tic
% for m = 1:M
%     Pi_for_moments = accumarray( grp2idx( num2str(Grouping_mat)), Pi_cond_mat(:,m), [], @mean) ;
%     G_for_moments = accumarray( grp2idx( num2str(Grouping_mat)), G_cond_mat(:,m), [], @mean) ;
% end
% toc

num_obs_mat = grpstats( repmat(C_mat,1,3) , {T_vec}, @sum);
Pi_for_moments = grpstats( Pi_cond_mat, {T_vec}, @sum) ./ num_obs_mat;
G_for_moments = grpstats( G_cond_mat, {T_vec}, @sum) ./ num_obs_mat;
Const_for_moments = grpstats( Const_cond_mat, {T_vec}, @sum) ./ num_obs_mat;

%Replace NaNs with 0s. This is the case when there are no firms that meet C
Pi_for_moments( isnan(Pi_for_moments)) = 0;
G_for_moments( isnan(G_for_moments)) = 0;
Const_for_moments( isnan(Const_for_moments)) = 0;

%Create a matrix with the values that will multiply theta_c
% For the first four moments, this is +1 if g>0 and -1 if g<0
% For the last two moents, this is 0

%Const_for_moments = sign(G_for_moments);
%Const_for_moments = repmat(sign(G_for_moments(:,1:6)),1,3);

%Set const to 0 for the last two types of moments
% Const_for_moments(:,5:6) = 0;
% Const_for_moments(:,11:12) = 0;
% Const_for_moments(:,17:18) = 0;

%Create a matrix that indicates whether each moment should be multiplied by
%lambda or not
%This is a matrix of size #grps x M where each column is all 1s if the mth moment
%relates to a case where J_tminus1= 1 and is 0 otherwise

Lambda_indicator_mat = repmat( moment_has_lambda , size(Pi_for_moments,1) , 3);


%Create a  fn that returns a (F*J) x M matrix where the (i,j)th entry is the jth moment for
%grp i

Interacted_Moments_mat_fn = @(theta_c, theta_g, lambda) (Pi_for_moments + ( lambda * Lambda_indicator_mat ...
    + (1 - Lambda_indicator_mat) ) .* ( theta_c * Const_for_moments + theta_g * G_for_moments));

Y = Pi_for_moments;
A_g = @(lambda)  (lambda * Lambda_indicator_mat  + (1 - Lambda_indicator_mat) ) .* G_for_moments ;
A_c = @(lambda)  ( lambda * Lambda_indicator_mat  + (1 - Lambda_indicator_mat) ) .* Const_for_moments ;

first_six_columns = @(mat) mat(:,1:6);

Moments_mat_fn = @(theta_c,theta_g, lambda) first_six_columns( Interacted_Moments_mat_fn(theta_c,theta_g,lambda));

Y_basic = first_six_columns(Y);
A_g_basic = @(lambda) first_six_columns(A_g(lambda));
A_c_basic = @(lambda) first_six_columns(A_c(lambda));
end
