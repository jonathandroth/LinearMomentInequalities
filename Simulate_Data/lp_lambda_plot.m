 ds_dir = strcat( data_output_dir, 'Interacted_Moments/');
 load(  strcat(ds_dir, 'lambda_identified_set') );
 load(  strcat(ds_dir, 'lambda_results') );
 
  
 plot(lambda_vec, [mean(conditional_test);...
                  mean(hybrid_test);...
                  mean(lf_test_original);...
                  mean(lf_test_modified)] ) 


identified_set_max = max( lambda_vec( lambda_identified_set == 1) );
identified_set_min = min( lambda_vec( lambda_identified_set == 1) );


if(~isempty(identified_set_max) )
    line( [identified_set_max; identified_set_max], [0;1], 'LineStyle', '--', 'Color',  'r');
else
    warning('Didnt find any lambdas in identfied set');
end

if(~isempty(identified_set_min) && identified_set_min ~= min(lambda_vec)  )
    line( [identified_set_min; identified_set_min], [0;1], 'LineStyle', '--', 'Color',  'r');
end

legend( 'Conditional', 'Hybrid', 'LF','LFN', 'Identified Set Bound', 'Location','eastoutside' );
ylabel('Rejection Probability');
xlabel('Lambda');
title('Rejection Probabilities for Lambda');
saveas( gcf, strcat(figures_output_dir,filename_graph ), 'epsc');

display( strcat(figures_output_dir,filename_graph ) )
