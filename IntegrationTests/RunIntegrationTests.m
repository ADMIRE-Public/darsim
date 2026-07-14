function RunIntegrationTests(logFileName)
%RUNINTEGRATIONTESTS Run all DARSim integration tests and write one log.
%
% Place this file in DARSim/IntegrationTests. Place the supplied helper
% files in DARSim/IntegrationTests/helper_scripts.

arguments
    logFileName string = "IntegrationTestLog.txt"
end

%% User-defined constraints
absoluteTolerance = 1.0e-10;
relativeTolerance = 1.0e-8;
deleteOutputAfterCompletion = true;

%% Locate required directories
integrationTestsDir = fileparts(mfilename('fullpath'));
darsimRoot = fileparts(integrationTestsDir);
srcDir = fullfile(darsimRoot, 'src');
fracDir = fullfile(srcDir, 'darsim_legacy', 'FracturesGenerator');
helper_scriptsDir = fullfile(integrationTestsDir, 'helper_scripts');

assert(isfolder(srcDir), 'Expected source directory not found: %s', srcDir);
assert(isfolder(helper_scriptsDir), 'Expected helper_scripts helper directory not found: %s', helper_scriptsDir);

originalDirectory = pwd;
addpath(helper_scriptsDir);
cleanupDirectory = onCleanup(@() cd(originalDirectory)); 
cleanupPath = onCleanup(@() localRemovePath(helper_scriptsDir));

logPath = fullfile(integrationTestsDir, logFileName);
logFileID = fopen(logPath, 'w');
assert(logFileID >= 0, 'Unable to open log file: %s', logPath);

suiteStart = tic;
writeLogHeader(logFileID, darsimRoot, srcDir, integrationTestsDir, ...
    absoluteTolerance, relativeTolerance, deleteOutputAfterCompletion);

%% Discover tests
info = discoverIntegrationTests(integrationTestsDir);
missingBenchmarks = findMissingBenchmarks(integrationTestsDir, info.existing);

fprintf(logFileID, 'Highest test number             : %d\n', info.highest);
fprintf(logFileID, 'Existing test directories       : %d\n', numel(info.existing));
fprintf(logFileID, 'Missing test directories        : %s\n', formatTestList(info.missing));
fprintf(logFileID, 'Missing benchmark directories   : %s\n\n', formatTestList(missingBenchmarks));

%% Pre-run cleanup
fprintf(logFileID, 'PRE-RUN CLEANUP\n%s\n', repmat('-', 1, 72));
deletedBefore = deleteOutputDirectories(integrationTestsDir, info.existing, logFileID);
fprintf(logFileID, 'Deleted %d pre-existing Output directories.\n\n', deletedBefore);

statuses = strings(info.highest, 1);
wallTimes = nan(info.highest, 1);

%% Run tests in numerical order
for testNumber = 1:info.highest
    testName = sprintf('Test%02d', testNumber);
    testDirectory = fullfile(integrationTestsDir, testName);

    fprintf(logFileID, '%s\n%s\n%s\n', repmat('=', 1, 72), testName, repmat('=', 1, 72));

    if ~isfolder(testDirectory)
        statuses(testNumber) = "skipped: missing directory";
        fprintf(logFileID, 'Status : %s\n\n', statuses(testNumber));
        continue;
    end

    if ~isfolder(fullfile(testDirectory, 'Output_benchmark'))
        statuses(testNumber) = "skipped: missing benchmark";
        fprintf(logFileID, 'Status : %s\n\n', statuses(testNumber));
        continue;
    end

    fprintf(logFileID, 'Status : running\n');
    testStart = tic;

    % DARSim may call fclose('all'), so do not keep the integration-test
    % log open while the simulator runs.
    fclose(logFileID);
    logFileID = -1;

    try
        % Check if fracture output needs to be generated
        fileId = fopen(append(testDirectory,'/Sim.txt'), 'r');
        fracLineID = 0;
        gridLineID = 0;
        while ~feof(fileId)
            line = fgetl(fileId);
            if fracLineID == 1
                testFRACTURED = line;
                fracLineID = 0;
            elseif gridLineID == 1
                testGRID = line;
                gridLineID = 0;
            end
            if contains(line, 'FRACTURED', 'IgnoreCase', false)
                fracLineID = 1;
            elseif contains(line, 'GRID_DISCRETIZATION', 'IgnoreCase', false)
                gridLineID = 1;
            end
        end

        % Create fracture output if required
        if strcmp(strtrim(testFRACTURED),'ON') && (strcmp(strtrim(testGRID),'CartesianGrid') || strcmp(strtrim(testGRID),'CornerPointGrid'))
            cd(fracDir)
            DARSimFracGen(testDirectory, 'FractureInput.txt');
        end
        
        % Run simulation
        cd(srcDir);
        DARSim(testDirectory, 'Sim.txt');
        wallTimes(testNumber) = toc(testStart);
    catch ME
        wallTimes(testNumber) = toc(testStart);
        cd(integrationTestsDir);
        logFileID = reopenLogForAppend(logPath);
        statuses(testNumber) = "failed: error";
        fprintf(logFileID, 'Status : %s\nWall time : %.6f s\n', statuses(testNumber), wallTimes(testNumber));
        writeExceptionToLog(logFileID, ME);
        fprintf(logFileID, '\n');
        continue;
    end

    cd(integrationTestsDir);
    logFileID = reopenLogForAppend(logPath);

    diaryPath = fullfile(testDirectory, 'Output', 'RunDiary.txt');

    diaryResult = parseRunDiary(diaryPath);

    if diaryResult.convergenceFailure
        statuses(testNumber) = "failed: convergence";
        fprintf(logFileID, 'Status : %s\nWall time : %.6f s\nDetails : %s\n\n', ...
            statuses(testNumber), wallTimes(testNumber), diaryResult.message);
        continue;
    elseif ~diaryResult.success
        statuses(testNumber) = "failed: error";
        fprintf(logFileID, 'Status : %s\nWall time : %.6f s\nDetails : %s\n\n', ...
            statuses(testNumber), wallTimes(testNumber), diaryResult.message);
        continue;
    end

    comparison = compareTestOutputs(testDirectory, testNumber, ...
        absoluteTolerance, relativeTolerance);

    if comparison.matches
        statuses(testNumber) = "pass";
        fprintf(logFileID, 'Status : %s\nWall time : %.6f s\nCompared VTK files : %d\n\n', ...
            statuses(testNumber), wallTimes(testNumber), comparison.comparedFileCount);
    else
        statuses(testNumber) = "failed: output mismatch";
        fprintf(logFileID, 'Status : %s\nWall time : %.6f s\n', ...
            statuses(testNumber), wallTimes(testNumber));
        writeComparisonFailureToLog(logFileID, comparison);
        fprintf(logFileID, '\n');
    end
end

%% Summary
summary = summarizeStatuses(statuses);
fprintf(logFileID, '%s\nSUMMARY\n%s\n', repmat('=', 1, 72), repmat('=', 1, 72));
fprintf(logFileID, 'Highest test number             : %d\n', info.highest);
fprintf(logFileID, 'Existing test directories       : %d\n', numel(info.existing));
fprintf(logFileID, 'Skipped: missing directory      : %d\n', summary.skippedMissingDirectory);
fprintf(logFileID, 'Skipped: missing benchmark      : %d\n', summary.skippedMissingBenchmark);
fprintf(logFileID, 'Passed                          : %d\n', summary.passed);
fprintf(logFileID, 'Failed: error                   : %d\n', summary.failedError);
fprintf(logFileID, 'Failed: convergence             : %d\n', summary.failedConvergence);
fprintf(logFileID, 'Failed: output mismatch         : %d\n', summary.failedOutputMismatch);
fprintf(logFileID, 'Total suite wall time           : %.6f s\n', toc(suiteStart));

if summary.failedError + summary.failedConvergence + summary.failedOutputMismatch == 0
    fprintf(logFileID, 'Overall status                  : PASS\n');
else
    fprintf(logFileID, 'Overall status                  : FAIL\n');
end

%% Post-run cleanup
fprintf(logFileID, '\nPOST-RUN CLEANUP\n%s\n', repmat('-', 1, 72));
if deleteOutputAfterCompletion
    deletedAfter = deleteOutputDirectories(integrationTestsDir, info.existing, logFileID);
    fprintf(logFileID, 'Deleted %d generated Output directories.\n', deletedAfter);
else
    fprintf(logFileID, 'Generated Output directories retained by user setting.\n');
end
fprintf(logFileID, '\nLog file : %s\n', logPath);
fclose(logFileID);
logFileID = -1;
end

function fileID = reopenLogForAppend(logPath)
fileID = fopen(logPath, 'a');
if fileID < 0
    error('DARSim integration tests:LogReopenFailed', ...
        'Unable to reopen integration-test log: %s', logPath);
end
end

function localRemovePath(folder)
try rmpath(folder); catch, end
end
