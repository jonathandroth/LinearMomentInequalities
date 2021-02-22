%filename_graph = 'Theta_g_Rejection_Probabilities';
if(use_zero_cutoff == 1)
    suffix = "_zerocutoff";
else
    suffix = '';
end

for filename_graph = {'Mean_Weight_Rejection_Probabilities','Theta_g_Rejection_Probabilities'}

filename_graph = filename_graph{1};

mat_05 = [];
mat_50 = [];
mat_95 = [];
mat_mean = [];


mat_05_coxandshi = [];
mat_50_coxandshi = [];
mat_95_coxandshi = [];
mat_mean_coxandshi = [];


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
excess_length_table_file = strcat(data_output_folder, filename_graph,'_', 'rows_for_excess_length_table', suffix);
load(excess_length_table_file)

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
    
filename_vec = [filename_vec ; filename_graph];
moment_type_cell{end+1} = moment_type{1}; 
    end
end

%%
num_params_column = [2;2;4;4;10;10];
num_moments_column = [6;14;14;38;38;110];

mat_05 = [ num_params_column, num_moments_column, mat_05];
mat_50 = [ num_params_column, num_moments_column, mat_50];
mat_95 = [ num_params_column, num_moments_column, mat_95];
mat_mean = [ num_params_column, num_moments_column, mat_mean];


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table','05',suffix,'.tex') ,'wt');
fprintf(fid,clean_latex(mat_05));
fclose(fid);


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table','50',suffix,'.tex') ,'wt');
fprintf(fid,clean_latex(mat_50));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table','95',suffix,'.tex') ,'wt');
fprintf(fid,clean_latex(mat_95));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table','mean',suffix,'.tex') ,'wt');
fprintf(fid,clean_latex(mat_mean));
fclose(fid);


%Print the table with cox and shi too

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table','05',suffix, '_coxandshi', '.tex') ,'wt');
fprintf(fid,clean_latex([mat_05, mat_05_coxandshi]));
fclose(fid);


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table','50',suffix, '_coxandshi', '.tex') ,'wt');
fprintf(fid,clean_latex([mat_50, mat_50_coxandshi]));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table','95',suffix, '_coxandshi', '.tex') ,'wt');
fprintf(fid,clean_latex([mat_95, mat_95_coxandshi]));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_excess_length_table','mean',suffix, '_coxandshi', '.tex') ,'wt');
fprintf(fid,clean_latex([mat_mean, mat_mean_coxandshi]));
fclose(fid);

end