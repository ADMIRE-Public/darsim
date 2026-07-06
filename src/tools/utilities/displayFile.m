function displayFile(fileName)
    % Open the file for reading
    fid = fopen(fileName, 'r');
    % Check if the file was opened successfully
    if fid == -1
        % Display an error message and exit the function
        fprintf('Error: could not open file %s\n', fileName);
        return;
    end
    % Read the file line by line and display each line
    while ~feof(fid)
        % Read a line from the file
        line = fgetl(fid);
        % Display the line
        disp(line);
    end
    % Close the file
    fclose(fid);
end