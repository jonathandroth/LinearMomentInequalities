%filename_graph = 'Theta_g_Rejection_Probabilities';
if(use_zero_cutoff == 1)
    suffix = "_zerocutoff";
else
    suffix = '';
end

for filename_graph = {'lambda_rejection_probabilities'}

filename_graph = filename_graph{1};

mat_05 = [];
mat_50 = [];
mat_95 = [];
mat_mean = [];

mat_05_with_infs = [];
mat_50_with_infs = [];
mat_95_with_infs = [];
mat_mean_with_infs = [];

mat_05_coxandshi = [];
mat_50_coxandshi = [];
mat_95_coxandshi = [];
mat_mean_coxandshi = [];

mat_05_with_infs_coxandshi = [];
mat_50_with_infs_coxandshi = [];
mat_95_with_infs_coxandshi = [];
mat_mean_with_infs_coxandshi = [];


filename_vec = [];
numthetha_vec = [];
moment_type_cell = {};


for numthetas = [1,3,9]
    for moment_type = {'Basic_Moments', 'Interacted_Moments'}
data_output_dir = strcat('../../Output/Conditional_FullMatrix/Data/',...
                             num2str(numthetas), 'Thetac', ifelse( numthetas==1, '','s'),...
                             '/', moment_type{1}, '/')

% data_output_dir = strcat('/Volumes/jonathanroth/Moment_Inequalities_Ariel/Output/Conditional_FullMatrix/Data/',...
%                              num2str(numthetas), 'Thetac', ifelse( numthetas==1, '','s'),...
%                              '/', moment_type{1}, '/')

%data_output_dir = '../../Output/Conditional_FullMatrix/Data/1Thetac/Basic_Moments/';

if( strcmp(filename_graph ,'Theta_g_Rejection_Probabilities') )
    data_output_dir = strcat(data_output_dir, 'Theta_g/');
end

data_output_folder = strcat( data_output_dir, dirname, 'Interacted_Moments/');
excess_length_table_lambda_file = strcat(data_output_folder, filename_graph,'_', 'rows_for_excess_length_table_lambda');
load(excess_length_table_lambda_file)

if(use_zero_cutoff == 0)
    mat_05 = [mat_05; row_05];
    mat_50 = [mat_50; row_50];
    mat_95 = [mat_95; row_95];
    mat_mean = [mat_mean; row_mean];
    
    mat_05_coxandshi = [mat_05_coxandshi; row_05_coxandshi];
    mat_50_coxandshi = [mat_50_coxandshi; row_50_coxandshi];
    mat_95_coxandshi = [mat_95_coxandshi; row_95_coxandshi];
    mat_mean_coxandshi = [mat_mean_coxandshi; row_mean_coxandshi];
else
    mat_05 = [mat_05; row_05_zerocutoff];
    mat_50 = [mat_50; row_50_zerocutoff];
    mat_95 = [mat_95; row_95_zerocutoff];
    mat_mean = [mat_mean; row_mean_zerocutoff];
    
    mat_05_coxandshi = [mat_05_coxandshi; row_05_zerocutoff_coxandshi];
    mat_50_coxandshi = [mat_50_coxandshi; row_50_zerocutoff_coxandshi];
    mat_95_coxandshi = [mat_95_coxandshi; row_95_zerocutoff_coxandshi];
    mat_mean_coxandshi = [mat_mean_coxandshi; row_mean_zerocutoff_coxandshi];
end

mat_05_with_infs = [mat_05_with_infs; row_05_with_infs];
mat_50_with_infs = [mat_50_with_infs; row_50_with_infs];
mat_95_with_infs = [mat_95_with_infs; row_95_with_infs];
mat_mean_with_infs = [mat_mean_with_infs; row_mean_with_infs];

mat_05_with_infs_coxandshi = [mat_05_with_infs_coxandshi; row_05_with_infs_coxandshi];
mat_50_with_infs_coxandshi = [mat_50_with_infs_coxandshi; row_50_with_infs_coxandshi];
mat_95_with_infs_coxandshi = [mat_95_with_infs_coxandshi; row_95_with_infs_coxandshi];
mat_mean_with_infs_coxandshi = [mat_mean_with_infs_coxandshi; row_mean_with_infs_coxandshi];


filename_vec = [filename_vec ; filename_graph];
moment_type_cell{end+1} = moment_type{1}; 
    end
end

%%

%Create a matrix that's one for entries that are affected by infinities
annotate_mat_05 = (Inf == mat_05_with_infs);
annotate_mat_50 = (Inf == mat_50_with_infs);
annotate_mat_95 = (Inf == mat_95_with_infs);
annotate_mat_mean = (Inf == mat_mean_with_infs);


annotate_mat_05_coxandshi = (Inf == mat_05_with_infs_coxandshi);
annotate_mat_50_coxandshi = (Inf == mat_50_with_infs_coxandshi);
annotate_mat_95_coxandshi = (Inf == mat_95_with_infs_coxandshi);
annotate_mat_mean_coxandshi = (Inf == mat_mean_with_infs_coxandshi);


mat_05 = [ [2;2;4;4;10;10]+1, [6;14;14;38;38;110], mat_05];
mat_50 = [ [2;2;4;4;10;10]+1, [6;14;14;38;38;110], mat_50];
mat_95 = [ [2;2;4;4;10;10]+1, [6;14;14;38;38;110], mat_95];
mat_mean = [ [2;2;4;4;10;10]+1, [6;14;14;38;38;110], mat_mean];

annotate_mat_05 = [ [0;0;0;0;0;0], [0;0;0;0;0;0], annotate_mat_05];
annotate_mat_50 = [ [0;0;0;0;0;0], [0;0;0;0;0;0], annotate_mat_50];
annotate_mat_95 = [ [0;0;0;0;0;0], [0;0;0;0;0;0], annotate_mat_95];
annotate_mat_mean = [ [0;0;0;0;0;0], [0;0;0;0;0;0], annotate_mat_mean];


annotationChar = '+';

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table_lambda','05',suffix,'.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_05, annotate_mat_05, annotationChar) ) );
fclose(fid);


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table_lambda','50',suffix,'.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_50, annotate_mat_50, annotationChar) ) );
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table_lambda','95',suffix,'.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_95, annotate_mat_95, annotationChar) ) );
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table_lambda','mean',suffix,'.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_mean, annotate_mat_mean, annotationChar) ) );

fclose(fid);



fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table_lambda','05',suffix, '_coxandshi', '.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_05_coxandshi, annotate_mat_05_coxandshi, annotationChar) ) );
fclose(fid);


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table_lambda','50',suffix, '_coxandshi', '.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_50_coxandshi, annotate_mat_50_coxandshi, annotationChar) ) );
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table_lambda','95',suffix, '_coxandshi', '.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_95_coxandshi, annotate_mat_95_coxandshi, annotationChar) ) );
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table_lambda','mean',suffix, '_coxandshi', '.tex') ,'wt');
fprintf( fid, clean_latex( mat2latexwithannotation( mat_mean_coxandshi, annotate_mat_mean_coxandshi, annotationChar) ) );
fclose(fid);



end