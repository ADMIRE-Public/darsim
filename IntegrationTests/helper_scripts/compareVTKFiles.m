function result=compareVTKFiles(newPath,benchPath,names,atol,rtol)
result=struct('matches',false,'messages',{{}},'variableResults',struct([]));
try new=readDARSimBinaryVTK(newPath,names); catch ME, result.messages{end+1}=sprintf('Cannot read new file: %s',ME.message); return; end
try bench=readDARSimBinaryVTK(benchPath,names); catch ME, result.messages{end+1}=sprintf('Cannot read benchmark: %s',ME.message); return; end
if new.cellCount~=bench.cellCount
    result.messages{end+1}=sprintf('CELL_DATA count differs: new=%d, benchmark=%d.',new.cellCount,bench.cellCount);
end
vr=repmat(struct('name','','matches',false,'message','','maxAbsoluteError',NaN,'maxRelativeError',NaN,'rmsError',NaN,'failingValueCount',0),numel(names),1);
for k=1:numel(names)
    name=names{k}; vr(k).name=name;
    hn=isfield(new.variables,name); hb=isfield(bench.variables,name);
    if ~hn||~hb
        if ~hn&&~hb, vr(k).message='Variable missing from both files.';
        elseif ~hn, vr(k).message='Variable missing from new file.';
        else, vr(k).message='Variable missing from benchmark file.'; end
        continue;
    end
    c=compareArrays(new.variables.(name),bench.variables.(name),atol,rtol);
    fields=fieldnames(c); for j=1:numel(fields), vr(k).(fields{j})=c.(fields{j}); end
end
result.variableResults=vr; result.matches=isempty(result.messages)&&all([vr.matches]);
end

function c=compareArrays(a,b,atol,rtol)
c=struct('matches',false,'message','','maxAbsoluteError',NaN,'maxRelativeError',NaN,'rmsError',NaN,'failingValueCount',0);
a=a(:); b=b(:);
if numel(a)~=numel(b), c.message=sprintf('Value count differs: new=%d, benchmark=%d.',numel(a),numel(b)); return; end
sameNonfinite=(isnan(a)&isnan(b))|(isinf(a)&isinf(b)&sign(a)==sign(b)); finite=isfinite(a)&isfinite(b); invalid=~(sameNonfinite|finite);
err=nan(size(a)); allowed=nan(size(a)); err(finite)=abs(a(finite)-b(finite)); allowed(finite)=atol+rtol.*abs(b(finite));
fail=invalid; fail(finite)=err(finite)>allowed(finite); c.failingValueCount=nnz(fail);
if any(finite)
    e=err(finite); c.maxAbsoluteError=max(e); c.rmsError=sqrt(mean(e.^2));
    c.maxRelativeError=max(e./max(abs(b(finite)),atol));
else, c.maxAbsoluteError=0; c.maxRelativeError=0; c.rmsError=0; end
c.matches=~any(fail);
if c.matches, c.message='Values match within tolerance.'; else, c.message=sprintf('%d of %d values failed.',c.failingValueCount,numel(a)); end
end
