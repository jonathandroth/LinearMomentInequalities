
%This function takes an input matrix and produces cleaned latex output using
% the latex function from symbolic math.
%It truncates all numbers to two decimal places and removes trailing zeros
%from integers
function txt = clean_latex( mat)

txt = latex(mat);
txt = regexprep(txt,'\.(\d\d)\d(\d)*\$','\.$1\$'); %convert to two decimal places
txt = regexprep(txt,'\.00\$','\$'); %convert two 0's after decimal to nothing
txt = regexprep(txt,'\\','\\\\'); %double backslashes
txt = strcat(txt,'\\\\');%add double backslash at the end

end