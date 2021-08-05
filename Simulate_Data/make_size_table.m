if(use_zero_cutoff == 1)
    suffix = "_zerocutoff";
else
    suffix = '';
end

for filename_graph = {'Mean_Weight_Rejection_Probabilities','Theta_g_Rejection_Probabilities'};

filename_graph = filename_graph{1};

mat_ub_sizes = [];
mat_lb_sizes = [];
mat_max_sizes = [];
mat_max_sizes_coxandshi = [];
mat_max_sizes_asandkms = [];

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
size_table_file = strcat(data_output_folder, filename_graph,'_', 'upper_and_lower_bound_size',suffix,'.mat');
load(size_table_file)

if(use_zero_cutoff == 0 )
    
    mat_ub_sizes = round([mat_ub_sizes; upper_bound_sizes],2);
    mat_lb_sizes = round([mat_lb_sizes; lower_bound_sizes],2);
    mat_max_sizes = round([mat_max_sizes; max_size_in_id_set],2);
    mat_max_sizes_coxandshi = round([mat_max_sizes_coxandshi; max_size_in_id_set_coxandshi],2);
    
    if(numthetas < 9)
      mat_max_sizes_asandkms = round([mat_max_sizes_asandkms; max_size_in_id_set_asandkms],2);
    end
else
    mat_ub_sizes = round([mat_ub_sizes; upper_bound_sizes_zerocutoff],2);
    mat_lb_sizes = round([mat_lb_sizes; lower_bound_sizes_zerocutoff],2);
    mat_max_sizes = round([mat_max_sizes; max_size_in_id_set_zerocutoff],2);
    mat_max_sizes_coxandshi = round([mat_max_sizes_coxandshi; max_size_in_id_set_zerocutoff_coxandshi],2);
    
    if(numthetas < 9)
      mat_max_sizes_asandkms = round([mat_max_sizes_asandkms; max_size_in_id_set_zerocutoff_asandkms],2);
    end
end


filename_vec = [filename_vec ; filename_graph];
moment_type_cell{end+1} = moment_type{1}; 
    end
end

%%

mat_ub_sizes = [ [2;2;4;4;10;10], [6;14;14;38;38;110], mat_ub_sizes];
mat_lb_sizes = [ [2;2;4;4;10;10], [6;14;14;38;38;110], mat_lb_sizes];

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_upper_bound_sizes',suffix,'.tex') ,'wt');
fprintf(fid,clean_latex(mat_ub_sizes));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_lower_bound_sizes',suffix,'.tex') ,'wt');
fprintf(fid,clean_latex(mat_lb_sizes));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_max_size_in_id_set', suffix, '.tex') ,'wt');
fprintf(fid,clean_latex(mat_max_sizes));
fclose(fid);

fid = fopen(strcat('../../Output/',...
                    filename_graph,'_max_size_in_id_set', suffix, '_coxandshi', '.tex') ,'wt');
fprintf(fid,clean_latex([mat_max_sizes, mat_max_sizes_coxandshi]));
fclose(fid);


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_max_size_in_id_set', suffix, '_asandkms', '.tex') ,'wt');
fprintf(fid,clean_latex([mat_max_sizes(1:4,:),mat_max_sizes_asandkms]));
fclose(fid);


fid = fopen(strcat('../../Output/',...
                    filename_graph,'_max_size_in_id_set', suffix, '_allcombined', '.tex') ,'wt');
fprintf(fid,clean_latex([mat_max_sizes, mat_max_sizes_coxandshi, [mat_max_sizes_asandkms; NaN, NaN; NaN, NaN]]));
fclose(fid);

end