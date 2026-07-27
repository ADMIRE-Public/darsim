function result = compareTestOutputs(testDirectory, testNumber, atol, rtol)
%COMPARETESTOUTPUTS Compare all expected reservoir and fracture VTK files.
newDir = fullfile(testDirectory, 'Output', 'VTK');
benchDir = fullfile(testDirectory, 'Output_benchmark', 'VTK');
result = struct('matches',false,'comparedFileCount',0,'generalMessages',{{}},'fileFailures',struct([]));
if ~isfolder(newDir), result.generalMessages{end+1}=sprintf('New VTK directory missing: %s',newDir); return; end
if ~isfolder(benchDir), result.generalMessages{end+1}=sprintf('Benchmark VTK directory missing: %s',benchDir); return; end
newNames=collectNames(newDir,testNumber); benchNames=collectNames(benchDir,testNumber);
if numel(newNames)~=numel(benchNames)
    result.generalMessages{end+1}=sprintf('VTK file count differs: new=%d, benchmark=%d.',numel(newNames),numel(benchNames));
end
missing=setdiff(benchNames,newNames); extra=setdiff(newNames,benchNames);
if ~isempty(missing), result.generalMessages{end+1}=sprintf('Missing new files: %s',strjoin(missing,', ')); end
if ~isempty(extra), result.generalMessages{end+1}=sprintf('Unexpected new files: %s',strjoin(extra,', ')); end
common=intersect(benchNames,newNames,'stable'); failures=struct([]);
for k=1:numel(common)
    c=compareVTKFiles(fullfile(newDir,common{k}),fullfile(benchDir,common{k}),{'P_1','P_2','S_1','S_2'},atol,rtol);
    result.comparedFileCount=result.comparedFileCount+1;
    if ~c.matches
        failures(end+1).fileName=common{k};
        failures(end).comparison=c;
    end
end
result.fileFailures=failures;
result.matches=isempty(result.generalMessages)&&isempty(failures);
end

function names=collectNames(folder,testNumber)
files = dir(fullfile(folder, '**', '*'));
files = files(~[files.isdir]);
allNames = {files.name};

idx = contains(allNames, "IntegrationTest_Fracture") | ...
      (contains(allNames, "IntegrationTest_Reservoir") & ...
      ~(contains(allNames, "Level")));

matchingFiles = files(idx);
names = {matchingFiles.name};
names=sort(unique(names));
end
