%% Set parameters for the lp confidence sets scripts
%To be called at the beginning of the run file


basic_inequalities_set_parameters   

diagonal = 0;

%working_dir = '/Volumes/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data';

%If don't have SLUM_CPUS_PER_TASK, then we are on laptop, set wd and
%numdatasets accordingly
if( isempty(getenv('SLURM_CPUS_PER_TASK')) )
   
   if(exist('remote_to_server') && remote_to_server == 1)
       working_dir = '/Volumes/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data/';
   else
       working_dir = '/Users/jonathanroth/Google Drive/Research Projects/Moment_Inequalities_Ariel/Code/Simulate_Data';
   end
   numdatasets = 2;
   onLaptop =1;

    
else
   working_dir = '/n/home12/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data';
   parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));
   numdatasets = 500;
   onLaptop =0 ;
end

cd( working_dir);
%Specify where the input data is (can be relative to the working_dir)
data_input_dir = '../../Output/Simulated_Data/';


mkdir(data_output_dir);
mkdir(figures_output_dir);

dirname = 'Calibrated_SigmaZeta/';


if( exist( 'F_group_cell_parameters') == 0)
    F_group_cell_parameters = F_group_cell_moments;
end


num_F_groups_moments = size(F_group_cell_moments,1);
num_F_groups_parameters = size(F_group_cell_parameters,1);


numsims_lp = 200;


if( exist( 'combine_theta_g_moments') == 0)
    combine_theta_g_moments = 0;
end

if( exist( 'use_basic_moments') == 0)
    use_basic_moments = 0;
end


if( exist( 'graphs_only') == 0 )
    graphs_only = 0;
end