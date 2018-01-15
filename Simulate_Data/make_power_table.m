%filename_graph = 'Theta_g_Rejection_Probabilities';
for filename_graph = {'Mean_Weight_Rejection_Probabilities','Theta_g_Rejection_Probabilities'};

filename_graph = filename_graph{1};

ub_mat_05 = [];
lb_mat_05 = [];
ub_mat_50 = [];
lb_mat_50 = [];
ub_mat_95 = [];
lb_mat_95 = [];

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
power_table_file = strcat(data_output_folder, filename_graph,'_', 'rows_for_power_table');
load(power_table_file)

ub_mat_05 = [ub_mat_05; ub_row_05];
lb_mat_05 = [lb_mat_05; lb_row_05];
ub_mat_50 = [ub_mat_50; ub_row_50];
lb_mat_50 = [lb_mat_50; lb_row_50];
ub_mat_95 = [ub_mat_95; ub_row_95];
lb_mat_95 = [lb_mat_95; lb_row_95];

filename_vec = [filename_vec ; filename_graph];
moment_type_cell{end+1} = moment_type{1}; 
    end
end

%%

mat05 = [ [2;2;4;4;10;10], [6;14;14;38;38;110],lb_mat_05,ub_mat_05];
mat50 = [ [2;2;4;4;10;10], [6;14;14;38;38;110],lb_mat_50,ub_mat_50];
mat95 = [ [2;2;4;4;10;10], [6;14;14;38;38;110],lb_mat_95,ub_mat_95];


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_power_table','05','.tex') ,'wt');
fprintf(fid,clean_latex(mat05));
fclose(fid);


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_power_table','50','.tex') ,'wt');
fprintf(fid,clean_latex(mat50));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_power_table','95','.tex') ,'wt');
fprintf(fid,clean_latex(mat95));
fclose(fid);

end