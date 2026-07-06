function result = findTwoLowestRowWise(matrix)
    % Step 1: Find the two lowest values in each row of the matrix
    lowestValues = zeros(size(matrix));
    for row = 1:size(matrix, 1)
        sortedRow = sort(matrix(row, :));
        lowestValues(row, :) = sortedRow(1:2);
    end

    % Step 2: Create a logical matrix of the same size
    result = false(size(matrix));

    % Step 3: Set true at the positions of the two lowest values in each row
    for row = 1:size(matrix, 1)
        result(row, :) = (matrix(row, :) == lowestValues(row, 1)) | (matrix(row, :) == lowestValues(row, 2));
    end
end