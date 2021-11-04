numgridpoints = 21;

%lambda_grid = linspace(0,1,numgridpoints);
theta_c_grid = linspace(80,180,numgridpoints);
theta_g_grid = linspace(-75,25, numgridpoints);

in_identified_set = ones( numgridpoints);

lambda_true = 0.386;
theta_c_true = 129.73;
theta_g_true = -21.38;

% theta_c_grid = theta_c_true;
% theta_g_grid = theta_g_true;

for ds = 1:50
%for ds = 1:500

load( strcat('../../Output/Simulated_Data/ds', num2str(ds), '.mat' ));

Pi_array = Pi_star_array;

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

G_gplus1_array = [ G_array(:,2:J,:), zeros(F,1,T)];
G_gminus1_array = [ zeros(F,1,T), G_array(:,1:(J-1),:)];

J_gp1_t_array = [ J_t_array(:,2:J,:), zeros(F,1,T)];
J_gm1_t_array = [ zeros(F,1,T), J_t_array(:,1:(J-1),:)];

J_gp1_tm1_array = [ J_tminus1_array(:,2:J,:), zeros(F,1,T)];
J_gm1_tm1_array = [ zeros(F,1,T), J_tminus1_array(:,1:(J-1),:)];

%Put into "long" format
Pi_vec = Pi_array(:);
Pi_gplus1_vec = Pi_gplus1_array(:);
Pi_gminus1_vec = Pi_gminus1_array(:);

% %Create a J_array for g+1
% J_gp1_t_array = [ J_t_array(:,(2:size(J_t_array,2)),:) , zeros(size(J_t_array,1) , 1, size(J_t_array,3) ) ];
% J_gm1_t_array = [ zeros(size(J_t_array,1) , 1, size(J_t_array,3) ) , J_t_array(:,(1:(size(J_t_array,2)-1) ),:) ];
% J_gp1_tm1_array = [ J_tminus1_array(:,(2:size(J_tminus1_array,2)),:) , zeros(size(J_tminus1_array,1) , 1, size(J_tminus1_array,3) ) ];
% J_gm1_tm1_array = [ zeros(size(J_tminus1_array,1) , 1, size(J_tminus1_array,3) ) , J_tminus1_array(:,(1:(size(J_tminus1_array,2)-1) ),:) ];


J_t_vec = J_t_array(:);
J_tminus1_vec =  J_tminus1_array(:);

J_gp1_t_vec = J_gp1_t_array(:);
J_gm1_t_vec = J_gm1_t_array(:);
J_gp1_tm1_vec = J_gp1_tm1_array(:);
J_gm1_tm1_vec = J_gm1_tm1_array(:);

G_vec = G_array(:);
G_gplus1_vec = G_gplus1_array(:);
G_gminus1_vec = G_gminus1_array(:);

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
C_mat(:,5) = (J_t_vec == 1 & J_tminus1_vec == 0 & G_vec ~= min(G_vec) & J_gm1_t_vec == 0 & J_gm1_tm1_vec == 0) ;
C_mat(:,6) = (J_t_vec == 1 & J_tminus1_vec == 0 & G_vec ~= max(G_vec) & J_gp1_t_vec == 0 & J_gp1_tm1_vec == 0);

moment_has_lambda = [1,1,0,0,0,0];%vector indicating whether the moment "needs a lambda"
moment_sign_on_fc = [1,-1,1,-1,1,1];
% Create a matrix where the i,j th entry is the value of delta-pi for observation
% i if C_ij == 1, and is 0 otherwise

%For the first four moments, deltapi is pi if J_t =1 ,and -pi if J_t = 0
%the first 4 moments. For the 5th moment, this is (Pi_vec -
%Pi_gminus1_vec); for the 6th, this is (Pi_vec - Pi_gplus1_vec)

Pi_cond_mat = [Pi_vec, -Pi_vec, Pi_vec, -Pi_vec,...
    Pi_vec - Pi_gminus1_vec, Pi_vec - Pi_gplus1_vec] .* C_mat;

% Create a matrix where the i,j th entry is the value that multiplies theta_g in the moment for observation
% i if C_ij == 1, and is 0 otherwise

%For moments 1:4, the value that multiplies theta_g is +g if the product is
%in, and -g if it is out. For the last two, it is -1 and +1
G_cond_mat = [-G_vec, G_vec, -G_vec, G_vec, G_gminus1_vec - G_vec , G_gplus1_vec - G_vec]  .* C_mat;


% %Avg the values of Pi and G for each (t,g) group. All of the moments will
% %be of the form: M = Pi_for_moments - constant1 - constant2 * G_for_moments
% Grouping_mat = [ T_vec, G_vec];

% tic
% for m = 1:M
%     Pi_for_moments = accumarray( grp2idx( num2str(Grouping_mat)), Pi_cond_mat(:,m), [], @mean) ;
%     G_for_moments = accumarray( grp2idx( num2str(Grouping_mat)), G_cond_mat(:,m), [], @mean) ;
% end
% toc

num_obs_mat = grpstats( C_mat, {T_vec, G_vec}, @sum);
Pi_for_moments = grpstats( Pi_cond_mat, {T_vec, G_vec}, @sum) ./ num_obs_mat;
G_for_moments = grpstats( G_cond_mat, {T_vec, G_vec}, @sum) ./ num_obs_mat;

%Replace NaNs with 0s. This is the case when there are no firms that meet C
Pi_for_moments( isnan(Pi_for_moments)) = 0;
G_for_moments( isnan(G_for_moments)) = 0;

%Create a matrix with the values that will multiply theta_c
% For the first four moments, this is +1 if g>0 and -1 if g<0
% For the last two moents, this is 0
Const_for_moments = sign(G_for_moments);
Const_for_moments(:,5:6) = 0;

%Create a matrix that indicates whether each moment should be multiplied by
%lambda or not
%This is a matrix of size #grps x M where each column is all 1s if the mth moment
%relates to a case where J_tminus1= 1 and is 0 otherwise

Lambda_indicator_mat = repmat( moment_has_lambda , size(Pi_for_moments,1) , 1);


%Create a  fn that returns a (F*J) x M matrix where the (i,j)th entry is the jth moment for
%grp i

%Also need to fix the grpstats so that they avg only over the things that
%are positive (i.e. the observations for which C > 0)
Moments_mat_fn = @(theta_c, theta_g, lambda) (Pi_for_moments + ( lambda * Lambda_indicator_mat ...
    + (1 - Lambda_indicator_mat) ) .* ( theta_c * Const_for_moments + theta_g * G_for_moments));

moment_fn = @(theta_c, theta_g) -Moments_mat_fn(theta_c, theta_g, lambda_true);
identified_set_ds = NaN(numgridpoints);

           
for theta_c_index = 1:size(theta_c_grid, 2);
    for theta_g_index = 1:size(theta_g_grid, 2);
        
        theta_c = theta_c_grid(theta_c_index);
        theta_g = theta_g_grid(theta_g_index);
        
        moments_mat = moment_fn(theta_c, theta_g);
        
        identified_set_ds(theta_c_index, theta_g_index) = min( moments_mat(:) <= 0 );
        
        
    end
end

in_identified_set = min( identified_set_ds, in_identified_set);

end
