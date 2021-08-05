function txt = mat2latex(mat) 
    
numRows = size(mat,1);
numCols = size(mat,2);

txt = '';

for(row = 1:numRows)
    for(col = 1:numCols)
        
        %Convert the relevant matrix entry to a string and add to txt,
        %bracketed by $ signs
        txt = [txt, '$', num2str( mat(row,col), '%.2f' ), '$' ];
        
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