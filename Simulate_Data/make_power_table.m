%filename_graph = 'Theta_g_Rejection_Probabilities';
filename_graph = 'Mean_Weight_Rejection_Probabilities';

ub_mat_05 = [];
lb_mat_05 = [];
ub_mat_50 = [];
lb_mat_50 = [];
ub_mat_95 = [];
lb_mat_95 = [];

filename_vec = [];
numthetha_vec = [];

for numthetas = [1,3,9]
data_output_dir = strcat('../../Output/Conditional_FullMatrix/Data/',...
                             num2str(numthetas), 'Thetac', ifelse( numthetas==1, '','s'),...
                             '/Basic_Moments/')

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

end

