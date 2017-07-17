    function [n m] = latex2(fmt,header,varargin)
%LATEX 
%  Creates a Latex table given fmt, a header, and a list of column vectors
%  to go in the table



%Set the number of columns to be the number of inputs-2
m = nargin-2;

mathstr = '$';

%Calculate the number of rows
n = 0;
for c = 1:m
   n = max(n, size(varargin{c},1)); 
end




%Print the top of the table
disp('\begin{table}[hbtp]')
disp('\begin{center}')

%Create m number of l's, for the tabular command
tabularcommand = '\begin{tabular}{';
for i = 1:m
   tabularcommand = strcat(tabularcommand, 'l '); 
end
tabularcommand =strcat(tabularcommand,'}');
disp(tabularcommand)





headertxt  = '';
%Loop through the elements of the header, and put &'s between them
for c = 1:m
   if c<m
        headertxt = strcat(headertxt, header{c},'&');
   else
       headertxt = strcat(headertxt, header{c});
   end
end
headertxt = strcat(headertxt,'\\');
disp(headertxt)
disp('\hline')

%Loop through the rows
for r = 1:n
    rowtext = '';
    %Loop through the columns
    for c = 1:m
       col = varargin{c};
       %If the column is shorter than the row number, make a blank entry
       if(size(col,1)<r)
           element = '';
           if c<m
           rowtext = strcat(rowtext, element,'&');
           else
           rowtext = strcat(rowtext, element);
           end
           %Otherwise, print the rth element of col
       else
           if ~iscell(col)
           element = num2str(col(r,:),fmt);
           else
           element = col{r};
           end
            %Add the element to the row text
            if c<m 
                if isnumeric(col(r,:))
                rowtext = strcat(rowtext, mathstr, element, mathstr,'&');
                else
                rowtext = strcat(rowtext, element, '&');   
                end
               
            else
                if isnumeric(col(r,:))
                rowtext = strcat(rowtext, mathstr, element, mathstr);
                else
                rowtext = strcat(rowtext, element);
                end
               
            
            end
       end
    end
     
      %At end of row, add a double slash
      rowtext = strcat(rowtext,'\\');
      disp(rowtext)  
end
disp( '\end{tabular}')
disp('\end{center}')
disp('\end{table}')
end