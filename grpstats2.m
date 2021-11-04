%This function is meant to emulate grpstats, but has better performance for
%large datasets

function results = grpstats2( data_mat , grp_vec) 

 % [~,~,pos] = unique(grp_mat,'rows');
% % Produce row, col subs 
 % [col,row] = meshgrid(1:size(data_mat,2),pos);
% % Accumulate 
 % result = [accumarray([row(:), col(:)], reshape(grp_mat,[],1),[],@mean)];


[~,~,pos] = unique(grp_vec,'rows');
% Produce row, col subs 
[col,row] = meshgrid(1:size(data_mat,2),pos);
% Accumulate 
results = [accumarray([row(:), col(:)], reshape(data_mat,[],1),[],@nanmean)];



end
