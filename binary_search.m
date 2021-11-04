function [x0 ,x0_vec, lb_vec, ub_vec, f_vec] = binary_search( f, xrange )

tol = 10^(-4);
maxiters = 10^5;

lb = xrange(1);
ub = xrange(2);

if( sign( f(lb) ) == sign( f(ub) )   )
    error( 'Sign at ub is same as sign at lb');
end


iters = 0;
fval = tol + 1;

x0_vec = NaN(1);
ub_vec = NaN(1);
lb_vec = NaN(1);
f_vec = NaN(1);

while ( abs( fval) > tol && iters< maxiters )
    
    x0 = 1/2 * (ub + lb);
    fval = f(x0);
 
    x0_vec = [x0_vec; x0];
    lb_vec = [lb_vec; lb];
    ub_vec = [ub_vec; ub];
    f_vec = [f_vec; fval];
    
       
    
    if( sign(fval) == sign( f(lb) ) )
        lb = x0;
    else
        ub = x0;
    end
    
    iters = iters + 1;
end

if( iters == maxiters )
    warning('Reached max iters without converging');

else
    display( strcat('Converged after  ', num2str(iters-1), ' iterations') ); 
end


end
