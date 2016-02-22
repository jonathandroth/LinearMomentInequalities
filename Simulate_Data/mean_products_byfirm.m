%This function takes an array that's assumed to be F x J x T
%It outputs an Fx1 vector indicating the mean number of products offered in a year
%for each firm, excluding the first BURNOUT years

function mean_by_firm = mean_products_byfirm( J_t_array, burnout )

%Remove the burnout years
J_t_wo_burnout = J_t_array( :,:,(burnout+1):(size(J_t_array,3)) );
%Sum along dimension two (products), then take mean along dim 3 (time)
mean_by_firm = mean(  sum( J_t_wo_burnout, 2) , 3);

end