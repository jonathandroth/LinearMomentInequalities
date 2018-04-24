function val = ifelse(condition, valIfTrue, valIfFalse)
    
    %Handle case where condition is a vector and valIfTrue of valIfFalse is
    %scalar
    if( length(valIfTrue) == 1 && length(condition)>1)
        valIfTrue = repmat( valIfTrue, size(condition) );
    end
       
    if( length(valIfFalse) == 1 && length(condition)>1)
        valIfFalse = repmat( valIfFalse, size(condition) );
    end
       
    
    val = valIfTrue;
    val( condition == false ) = valIfFalse(condition == false);


end


