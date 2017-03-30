panel = reshape( (1:27), 3,3,3);



z = [1;1;1;2;2;3];
x = [1;1;2;2;2;3];

%% This is what you want to use:
subset = panel(sub2ind( size(panel), repmat(x, 1, size(panel,2)), repmat( (1:size(panel,2)), size(x,1), 1) ,repmat(z, 1, size(panel,2)) ))


subset = panel(x,:,z)


panel([3; 2],2,[1;2])


mat = [1,2;3,4]

mat( 1:2, 1:2)

panel( [1 2 3; 4 5 6])
