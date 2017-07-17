
%This function takes values x1, x2, y1, y2, y_target where either
% y1<= y_target <= y2 or y2<=y_target<=y1, and finds the value x such that
% f(x) = y_target using a linear interpolation between (x1,y1) and (x2,y2)
%If y1=y2=y_target, returns x1
function x_target = x_interpolated( x1, x2, y1, y2, y_target)

if(y1 == y2 && y1 == y_target)
    x_target = x1;
    return;
else
    
t = (y_target - y2) ./ (y1 - y2);

if( ~ ( 0<=t && t<=1) )
    error( 'y_target not between 0 and 1');
end

x_target = t*x1 + (1-t) * x2;



end


end