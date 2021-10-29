
%This function takes an input matrix and produces cleaned latex output using
% the latex function from symbolic math.
%It truncates all numbers to two decimal places and removes trailing zeros
%from integers

%Alternatively, the input can be text, in which case it just cleans it

function txt = clean_latex( mat)

if(isnumeric(mat) )    
    %txt = latex(mat);
    txt = mat2latex(mat);
else
    txt = mat;
end

txt = regexprep(txt,'\.(\d\d)\d(\d)*\$','\.$1\$'); %convert to two decimal places
txt = regexprep(txt,'([1-9])\.00\$','$1\$'); %convert two 0's after decimal to nothing, unless 0.0
txt = regexprep(txt,'([1-9]0)\.00\$','$1\$'); %this fixes numbers like 10.00 to 10
txt = regexprep(txt,'\\','\\\\'); %double backslashes
txt = regexprep(txt, 'NaN', ' '); %NaNs to blanks
txt = strcat(txt,'\\\\');%add double backslash at the end

end