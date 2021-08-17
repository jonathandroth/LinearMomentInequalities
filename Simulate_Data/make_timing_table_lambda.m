
for filename_graph = {'lambda_rejection_probabilities'}


filename_graph = filename_graph{1};

mat_timing_allcombined = [];

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
timing_table_file = strcat(data_output_folder, filename_graph,'_', 'timing_row_allcombined','.mat');
load(timing_table_file)
    
mat_timing_allcombined = [mat_timing_allcombined; timing_row_allcombined];

filename_vec = [filename_vec ; filename_graph];
moment_type_cell{end+1} = moment_type{1}; 
    end
end

%%

%add row numbers and convert to minutes
mat_timing_allcombined = [ [3;3;5;5;11;11], [6;14;14;38;38;110], mat_timing_allcombined/60];


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_mat_timing_allcombined','.tex') ,'wt');
fprintf(fid,clean_latex(mat_timing_allcombined));
fclose(fid);


end