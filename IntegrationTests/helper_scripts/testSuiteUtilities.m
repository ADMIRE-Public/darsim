function varargout = testSuiteUtilities(action, varargin)
%TESTSUITEUTILITIES Dispatcher for small integration-test utilities.
switch action
    case 'discover', varargout{1} = discover(varargin{:});
    case 'missingBenchmarks', varargout{1} = missingBenchmarks(varargin{:});
    case 'deleteOutputs', varargout{1} = deleteOutputs(varargin{:});
    case 'formatList', varargout{1} = formatList(varargin{:});
    case 'header', writeHeader(varargin{:});
    case 'parseDiary', varargout{1} = parseDiary(varargin{:});
    case 'exception', writeException(varargin{:});
    case 'summary', varargout{1} = summary(varargin{:});
    otherwise, error('Unknown utility action: %s', action);
end
end

function info = discover(root)
entries = dir(fullfile(root, 'Test*'));
entries = entries([entries.isdir]);
numbers = [];
for k = 1:numel(entries)
    token = regexp(entries(k).name, '^Test(\d+)$', 'tokens', 'once');
    if ~isempty(token), numbers(end+1) = str2double(token{1}); end %#ok<AGROW>
end

if isempty(numbers)
    info = struct('existing', [], 'highest', 0, 'missing', []);
else
    info = struct('existing', numbers, 'highest', max(numbers), ...
        'missing', setdiff(1:max(numbers), numbers));
end
end

function missing = missingBenchmarks(root, numbers)
missing = [];
for n = numbers
    if ~isfolder(fullfile(root, sprintf('Test%02d', n), 'Output_benchmark'))
        missing(end+1) = n; %#ok<AGROW>
    end
end
end

function count = deleteOutputs(root, numbers, logID)
count = 0;
for n = numbers
    folder = fullfile(root, sprintf('Test%02d', n), 'Output');
    fracFile1 = fullfile(root, sprintf('Test%02d', n), 'Fracture_Output.txt');
    fracFile2 = fullfile(root, sprintf('Test%02d', n), 'fracture_process.mat');
    fracFolder = fullfile(root, sprintf('Test%02d', n), 'Output_FracGen');
    CPGFile = fullfile(root, sprintf('Test%02d',n), 'CornerPointGridData.mat');
    if isfolder(folder)
        try
            rmdir(folder, 's'); count = count + 1;
            fprintf(logID, 'Deleted : %s\n', folder);
        catch ME
            fprintf(logID, 'WARNING: Could not delete %s\nReason: %s\n', folder, ME.message);
        end
    end
    if isfile(fracFile1)
        try
            delete(fracFile1);
        catch ME
            fprintf(logID, 'WARNING: Could not delete %s\nReason: %s\n', fracFile1, ME.message);
        end
    end
    if isfile(fracFile2)
        try
            delete(fracFile2);
        catch ME
            fprintf(logID, 'WARNING: Could not delete %s\nReason: %s\n', fracFile2, ME.message);
        end
    end
    if isfolder(fracFolder)
        try
            rmdir(fracFolder, 's');
            fprintf(logID, 'Deleted : %s\n', fracFolder);
        catch ME
            fprintf(logID, 'WARNING: Could not delete %s\nReason: %s\n', fracFolder, ME.message);
        end
    end
    if isfile(CPGFile)
        try
            delete(CPGFile);
        catch ME
            fprintf(logID, 'WARNING: Could not delete %s\nReason: %s\n', CPGFile, ME.message);
        end
    end
end
end

function text = formatList(numbers)
if isempty(numbers), text = 'none'; return; end
parts = arrayfun(@(n) sprintf('Test%d', n), numbers, 'UniformOutput', false);
text = strjoin(parts, ', ');
end

function writeHeader(fileID, root, src, tests, atol, rtol, deleteAfter)
fprintf(fileID, '%s\nDARSim Integration Test Report\n%s\n', repmat('=',1,72), repmat('=',1,72));
fprintf(fileID, 'Date and time                   : %s\n', char(datetime('now')));
fprintf(fileID, 'MATLAB version                  : %s\n', version);
fprintf(fileID, 'DARSim root                     : %s\n', root);
fprintf(fileID, 'Source directory                : %s\n', src);
fprintf(fileID, 'Integration test directory      : %s\n', tests);
fprintf(fileID, 'Absolute tolerance              : %.16g\n', atol);
fprintf(fileID, 'Relative tolerance              : %.16g\n', rtol);
fprintf(fileID, 'Delete outputs after completion : %s\n\n', string(deleteAfter));
end

function result = parseDiary(path)
result = struct('success', false, 'convergenceFailure', false, 'message', '');
if ~isfile(path), result.message = sprintf('RunDiary.txt was not found: %s', path); return; end
try text = fileread(path); catch ME, result.message = ME.message; return; end
lines = splitlines(string(strrep(strrep(text, sprintf('\r\n'), sprintf('\n')), sprintf('\r'), sprintf('\n'))));
while ~isempty(lines) && strlength(strtrim(lines(end))) == 0, lines(end) = []; end
if isempty(lines), result.message = 'RunDiary.txt is empty.'; return; end
if strtrim(lines(end)) == "Simulation successfully finished!"
    result.success = true; return;
end
match = regexp(char(text), 'The solution has not converged at timestep[^\r\n]*', 'match', 'once');
if ~isempty(match)
    result.convergenceFailure = true; result.message = strtrim(match);
else
    result.message = sprintf('Final nonblank diary line was: "%s"', strtrim(lines(end)));
end
end

function writeException(fileID, ME)
fprintf(fileID, 'Error identifier : %s\nError message    : %s\n', ME.identifier, ME.message);
if ~isempty(ME.stack)
    fprintf(fileID, 'Stack trace:\n');
    for k=1:numel(ME.stack), fprintf(fileID, '  %s, line %d\n', ME.stack(k).name, ME.stack(k).line); end
end
end

function out = summary(values)
out = struct('passed', nnz(values=="pass"), ...
    'skippedMissingDirectory', nnz(values=="skipped: missing directory"), ...
    'skippedMissingBenchmark', nnz(values=="skipped: missing benchmark"), ...
    'failedError', nnz(values=="failed: error"), ...
    'failedConvergence', nnz(values=="failed: convergence"), ...
    'failedOutputMismatch', nnz(values=="failed: output mismatch"));
end