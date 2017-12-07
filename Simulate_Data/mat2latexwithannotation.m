function txt = mat2latexwithannotation(mat,annotationMat, annotationTxt) 
    
numRows = size(mat,1);
numCols = size(mat,2);

txt = '';

for(row = 1:numRows)
    for(col = 1:numCols)
        
        dataEntryString = ['$', num2str( mat(row,col) ), '$'];
        
        if( annotationMat(row,col) == 1)
            dataEntryString = [dataEntryString, '^{', annotationTxt, '}'];
        end
                
        %Convert the relevant matrix entry to a string and add to txt,
        %bracketed by $ signs
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