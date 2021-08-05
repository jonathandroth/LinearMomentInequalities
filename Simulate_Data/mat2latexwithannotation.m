function txt = mat2latexwithannotation(mat,annotationMat, annotationTxt) 
    
numRows = size(mat,1);
numCols = size(mat,2);

txt = '';

for(row = 1:numRows)
    for(col = 1:numCols)
        
        %Load the relevant entry in mat, round to two decimals, convert to string
        dataEntryString = num2str( num2str( mat(row,col), '%.2f' ) );
        
        dataEntryString = regexprep(dataEntryString,'([1-9])\.00','$1'); %convert two 0's after decimal to nothing, unless 0.0
        dataEntryString = regexprep(dataEntryString,'([1-9]0)\.00','$1'); %this fixes numbers like 10.00 to 10
        
        %If we should annotate, add the annotation as a superscript
        if( annotationMat(row,col) == 1)
            dataEntryString = [dataEntryString, '\\rlap{$^{', annotationTxt, '}$}'];
        end
        
     
        %Convert the relevant matrix entry to a string and add to txt
        
        %txt = [txt, ['$', dataEntryString, '$']];
        txt = [txt, dataEntryString];
        
        %If not last column, add spaces and alignment sign
        if( col ~= numCols)
            txt = [txt, ' & '];
        end
        
        %If last column, add ' \\ \n' if we're not in the last row
        if( col == numCols && row ~= numRows)
            txt = [txt, ' \\\\ \n'];
        end
        
    end
    

end
    txt = sprintf(txt);
end