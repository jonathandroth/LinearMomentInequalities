%% Set parameters for the lp confidence sets scripts
%To be called at the beginning of the run file


basic_inequalities_set_parameters   

diagonal = 0;

%working_dir = '/Volumes/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data';
%working_dir = '/Users/jonathanroth/Google Drive/Research Projects/Moment_Inequalities_Ariel/Code/Simulate_Data';
working_dir = '/n/home12/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data';
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

cd( working_dir);
%Specify where the input data is (can be relative to the working_dir)
data_input_dir = '../../Output/Simulated_Data/';


mkdir(data_output_dir);
mkdir(figures_output_dir);

dirname = 'Calibrated_SigmaZeta/';



num_F_groups = size(F_group_cell,1);

numsims_lp = 200;


if( exist( 'combine_theta_g_moments') == 0)
    combine_theta_g_moments = 0;
end

if( exist( 'use_basic_moments') == 0)
    use_basic_moments = 0;
end