function val = ifelse(condition, valIfTrue, valIfFalse)
    if condition
        val = valIfTrue;
    else
        val = valIfFalse;
    end
end