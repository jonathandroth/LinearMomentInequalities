
%Specify the working directory

%working_dir = '/Users/jonathanroth/Google Drive/Research Projects/Moment_Inequalities_Ariel/Code/Simulate_Data';
working_dir = '/n/home12/jonathanroth/Moment_Inequalities_Ariel/Code/Simulate_Data';
%parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

cd( working_dir);
%Specify where the input data is (can be relative to the working_dir)
data_input_dir = '../../Output/Simulated_Data/';

%Specify where the output should go (can be relative to the working  dir)
data_output_dir = '../../Output/Conditional_FullMatrix/Data/';
figures_output_dir = '../../Figures/Conditional_FullMatrix/';

%Specify the subdirectory names within the data folder (e.g. if have different DGPs). 
dirnames = { 'Calibrated_SigmaZeta/'};


simulate_data_main
