function result = cellReduce( c, fun)

numDimensions = length(size(c{1}));

concatenatedCell = cat(numDimensions +1, c{:} );

result = fun(concatenatedCell, numDimensions + 1);

end