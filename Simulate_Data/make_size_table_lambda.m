if(use_zero_cutoff == 1)
    suffix = "_zerocutoff";
else
    suffix = '';
end

for filename_graph = {'lambda_rejection_probabilities'}


filename_graph = filename_graph{1};

mat_ub_sizes = [];
mat_lb_sizes = [];

filename_vec = [];
numthetha_vec = [];
moment_type_cell = {};


for numthetas = [1,3,9]
    for moment_type = {'Basic_Moments', 'Interacted_Moments'}
%data_output_dir = strcat('../../Output/Conditional_FullMatrix/Data/',...
%                             num2str(numthetas), 'Thetac', ifelse( numthetas==1, '','s'),...
%                             '/', moment_type{1}, '/')

data_output_dir = strcat('../../Output/Conditional_FullMatrix/Data/',...
                             num2str(numthetas), 'Thetac', ifelse( numthetas==1, '','s'),...
                             '/', moment_type{1}, '/')

%data_output_dir = '../../Output/Conditional_FullMatrix/Data/1Thetac/Basic_Moments/';

if( strcmp(filename_graph ,'Theta_g_Rejection_Probabilities') )
    data_output_dir = strcat(data_output_dir, 'Theta_g/');
end

data_output_folder = strcat( data_output_dir, dirname, 'Interacted_Moments/');
size_table_file = strcat(data_output_folder, filename_graph,'_', 'upper_bound_size_lambda',suffix,'.mat');
load(size_table_file)

if(use_zero_cutoff == 0 )
    
    mat_ub_sizes = round([mat_ub_sizes; upper_bound_sizes],2);
    %mat_lb_sizes = round([mat_lb_sizes; lower_bound_sizes],2);
else
    mat_ub_sizes = round([mat_ub_sizes; upper_bound_sizes_zerocutoff],2);
    %mat_lb_sizes = round([mat_lb_sizes; lower_bound_sizes_zerocutoff],2);
end

filename_vec = [filename_vec ; filename_graph];
moment_type_cell{end+1} = moment_type{1}; 
    end
end

%%

mat_ub_sizes = [ [2;2;4;4;10;10], [6;14;14;38;38;110], mat_ub_sizes];
%mat_lb_sizes = [ [2;2;4;4;10;10], [6;14;14;38;38;110], mat_lb_sizes];

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_upper_bound_sizes_lambda',suffix,'.tex') ,'wt');
fprintf(fid,clean_latex(mat_ub_sizes));
fclose(fid);

% fid = fopen(strcat('../../Output/',...
%                     filename_graph,'_lower_bound_sizes',suffix,'.tex') ,'wt');
% fprintf(fid,clean_latex(mat_lb_sizes));
% fclose(fid);


end