function [heightMatrix] = readHeightMatrix(fileName)
    heightMatrix = dlmread(fileName, ' ');
    
    % Matlab insert one extra column (become 513x514 instead of 513x513).
    % Following line deletes that last column
    heightMatrix(:,size(heightMatrix,2)) = [];

    fprintf('Succesfully read file: %s\n', fileName);
    fprintf('Resolution of heightmap = %d x %d\n', size(heightMatrix,1), size(heightMatrix,2));
end