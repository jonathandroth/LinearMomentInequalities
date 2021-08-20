%% Specify parameters
lp_set_parameters



%% Import data and calculate moments and covariance matrices and accept/reject over grid
if( graphs_only == 0 )
    
    %If on laptop (doing timing), expand grid to 1000
    if(onLaptop == 1)
        lambda_vec = linspace(min(lambda_vec),max(lambda_vec),1000)';
    end
lp_lambda_moments_and_covariances    

%% Compute identified set

tic;

display('Began identified set calc');

lp_lambda_identified_set
    
display('Finished identified set calc');

toc;

end

 %% Plot

lp_lambda_plot