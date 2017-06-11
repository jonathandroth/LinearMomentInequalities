
cd('/Volumes/jonathanroth/Moments_Inequalities_Ariel/Figures/Conditional_FullMatrix')

figures_dir = '~/Dropbox/Linear Moment Inequalities/Figures/'

rejection_grids_dir = strcat(figures_dir, 'Rejection_Grids')
mkdir(rejection_grids_dir)

lp_figures_dir = strcat(figures_dir, 'LP_figures')
mkdir(lp_figures_dir)

copyfile( 'Calibrated_SigmaZeta/*' , rejection_grids_dir)

copyfile( 'LP_figures/' , lp_figures_dir)
