%filename_graph = 'Theta_g_Rejection_Probabilities';
for filename_graph = {'Mean_Weight_Rejection_Probabilities','Theta_g_Rejection_Probabilities'};

filename_graph = filename_graph{1};

mat_05 = [];
mat_50 = [];
mat_95 = [];
mat_mean = [];

mat_05_with_infs = [];
mat_50_with_infs = [];
mat_95_with_infs = [];
mat_mean_with_infs = [];

filename_vec = [];
numthetha_vec = [];
moment_type_cell = {};


for numthetas = [1,3,9]
    for moment_type = {'Basic_Moments', 'Interacted_Moments'}
%data_output_dir = strcat('../../Output/Conditional_FullMatrix/Data/',...
%                             num2str(numthetas), 'Thetac', ifelse( numthetas==1, '','s'),...
%                             '/', moment_type{1}, '/')

data_output_dir = strcat('/Volumes/jonathanroth/Moment_Inequalities_Ariel/Output/Conditional_FullMatrix/Data/',...
                             num2str(numthetas), 'Thetac', ifelse( numthetas==1, '','s'),...
                             '/', moment_type{1}, '/')

%data_output_dir = '../../Output/Conditional_FullMatrix/Data/1Thetac/Basic_Moments/';

if( strcmp(filename_graph ,'Theta_g_Rejection_Probabilities') )
    data_output_dir = strcat(data_output_dir, 'Theta_g/');
end

data_output_folder = strcat( data_output_dir, dirname, 'Interacted_Moments/');
excess_length_table_lambda_file = strcat(data_output_folder, filename_graph,'_', 'rows_for_excess_length_table_lambda');
load(excess_length_table_lambda_file)

mat_05 = [mat_05; row_05];
mat_50 = [mat_50; row_50];
mat_95 = [mat_95; row_95];
mat_mean = [mat_mean; row_mean];

mat_05_with_infs = [mat_05_with_infs; row_05_with_infs];
mat_50_with_infs = [mat_50_with_infs; row_50_with_infs];
mat_95_with_infs = [mat_95_with_infs; row_95_with_infs];
mat_mean_with_infs = [mat_mean_with_infs; row_mean_with_infs];

filename_vec = [filename_vec ; filename_graph];
moment_type_cell{end+1} = moment_type{1}; 
    end
end

%%

%Create a matrix that's one for entries that are affected by infinities
annonate_mat_05 = (mat_05 == mat_05_with_infs);
annonate_mat_50 = (mat_50 == mat_50_with_infs);
annonate_mat_95 = (mat_95 == mat_95_with_infs);
annonate_mat_mean = (mat_mean == mat_mean_with_infs);


mat_05 = [ [2;2;4;4;10;10]+1, [6;14;14;38;38;110], mat_05];
mat_50 = [ [2;2;4;4;10;10]+1, [6;14;14;38;38;110], mat_50];
mat_95 = [ [2;2;4;4;10;10]+1, [6;14;14;38;38;110], mat_95];
mat_mean = [ [2;2;4;4;10;10]+1, [6;14;14;38;38;110], mat_mean];

annotate_mat_05 = [ [0;0;0;0;0;0], [0;0;0;0;0;0], annotate_mat_05];
annotate_mat_50 = [ [0;0;0;0;0;0], [0;0;0;0;0;0], annotate_mat_50];
annotate_mat_95 = [ [0;0;0;0;0;0], [0;0;0;0;0;0], annotate_mat_95];
annotate_mat_mean = [ [0;0;0;0;0;0], [0;0;0;0;0;0], annotate_mat_mean];

annotationChar = '+';

fid = fopen(strcat('/Volumes/jonathanroth/Moment_Inequalities_Ariel/Output/',...
                    filename_graph,'_excess_length_table_lambda','05','.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_05, annotate_mat_05, annotationChar) ) );
fclose(fid);


fid = fopen(strcat('/Volumes/jonathanroth/Moment_Inequalities_Ariel/Output/',...
                    filename_graph,'_excess_length_table_lambda','50','.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_50, annotate_mat_50, annotationChar) ) );
fclose(fid);

fid = fopen(strcat('/Volumes/jonathanroth/Moment_Inequalities_Ariel/Output/',...
                    filename_graph,'_excess_length_table_lambda','95','.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_95, annotate_mat_95, annotationChar) ) );
fclose(fid);

fid = fopen(strcat('/Volumes/jonathanroth/Moment_Inequalities_Ariel/Output/',...
                    filename_graph,'_excess_length_table_lambda','mean','.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_mean, annotate_mat_mean, annotationChar) ) );

fclose(fid);

end