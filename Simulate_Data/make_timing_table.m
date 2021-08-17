
for filename_graph = {'Mean_Weight_Rejection_Probabilities','Theta_g_Rejection_Probabilities'}

filename_graph = filename_graph{1};

timing_mat = [];
timing_mat_coxandshi = [];
timing_mat_asandkms = [];

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
timing_table_file = strcat(data_output_folder, filename_graph,'_', 'timing_row');
load(timing_table_file)

timing_mat = [timing_mat; timing_row];
timing_mat_coxandshi = [timing_mat_coxandshi; timing_row_coxandshi];

if(numthetas < 9)
timing_mat_asandkms = [timing_mat_asandkms; timing_row_asandkms];
end

    
filename_vec = [filename_vec ; filename_graph];
moment_type_cell{end+1} = moment_type{1}; 
    end
end

%%
num_params_column = [2;2;4;4;10;10];
num_moments_column = [6;14;14;38;38;110];

%Convert to minutes
timing_mat = timing_mat / 60;
timing_mat_coxandshi = timing_mat_coxandshi / 60;
timing_mat_asandkms = timing_mat_asandkms / 60;


timing_mat = [ num_params_column, num_moments_column, timing_mat];

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_timing_table','.tex') ,'wt');
fprintf(fid,clean_latex(timing_mat));
fclose(fid);


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_timing_table_coxandshi','.tex') ,'wt');
fprintf(fid,clean_latex(timing_mat_coxandshi));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_timing_table_asandkms','.tex') ,'wt');
fprintf(fid,clean_latex(timing_mat_asandkms));
fclose(fid);


timing_mat_allcombined = [timing_mat, timing_mat_coxandshi, [timing_mat_asandkms; NaN(2,2)] ];

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_timing_table_allcombined','.tex') ,'wt');
fprintf(fid,clean_latex(timing_mat_allcombined));
fclose(fid);
end