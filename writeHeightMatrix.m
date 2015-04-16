function [] = writeHeightMatrix(heightMatrix, fileName)
    dlmwrite(fileName, heightMatrix, 'delimiter', ' ', 'precision', 8, 'newline', 'pc');   
    fprintf('Succesfully written to file\n');
end

