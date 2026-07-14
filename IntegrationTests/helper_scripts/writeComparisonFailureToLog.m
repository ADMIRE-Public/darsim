function writeComparisonFailureToLog(fileID,c)
for k=1:numel(c.generalMessages), fprintf(fileID,'  - %s\n',c.generalMessages{k}); end
for f=1:numel(c.fileFailures)
    x=c.fileFailures(f); fprintf(fileID,'VTK file : %s\n',x.fileName);
    for k=1:numel(x.comparison.messages), fprintf(fileID,'  - %s\n',x.comparison.messages{k}); end
    for k=1:numel(x.comparison.variableResults)
        v=x.comparison.variableResults(k); if v.matches, continue; end
        fprintf(fileID,'  Variable : %s\n    Details                : %s\n',v.name,v.message);
        if ~isnan(v.maxAbsoluteError), fprintf(fileID,'    Maximum absolute error : %.16g\n',v.maxAbsoluteError); end
        if ~isnan(v.maxRelativeError), fprintf(fileID,'    Maximum relative error : %.16g\n',v.maxRelativeError); end
        if ~isnan(v.rmsError), fprintf(fileID,'    RMS error              : %.16g\n',v.rmsError); end
        fprintf(fileID,'    Failing values         : %d\n',v.failingValueCount);
    end
end
end
