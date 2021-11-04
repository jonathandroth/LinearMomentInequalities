%This function returns x mod y, except instead of returning zeros for
%elements of x that are multiples of y, it returns y for those elements

% E.g mod( [1;2;3] , 3) = [1;2;0], but mod2([1;2;3],3) = [1;2;3]

function result = mod2(x,y)

result = mod(x,y);
result( result == 0) = y;

end
