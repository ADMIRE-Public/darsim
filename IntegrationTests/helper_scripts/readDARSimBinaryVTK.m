function output=readDARSimBinaryVTK(filePath,requested)
%READDARSIMBINARYVTK Read selected CELL_DATA scalars from DARSim binary VTK.
fid=fopen(filePath,'r','ieee-be');
if fid<0, error('Unable to open %s',filePath); end
clean=onCleanup(@() fclose(fid)); %#ok<NASGU>
line=nonemptyLine(fid);
if ~startsWith(strtrim(line),'# vtk DataFile Version'), error('Invalid legacy VTK header.'); end
fgetl(fid); format=nonemptyLine(fid);
if ~strcmpi(strtrim(format),'BINARY'), error('Only binary legacy VTK is supported.'); end
dataset=nonemptyLine(fid);
if isempty(regexp(strtrim(dataset),'^DATASET\s+\S+','once')), error('Missing DATASET declaration.'); end
cellCount=[]; vars=struct();
while ~feof(fid)
    line=fgetl(fid); if ~ischar(line), break; end
    line=strtrim(line); if isempty(line), continue; end
    t=strsplit(line); key=upper(t{1});
    switch key
        case {'DIMENSIONS'}
        case {'X_COORDINATES','Y_COORDINATES','Z_COORDINATES'}
            skipBlock(fid,toCount(t{2},line),vtkPrecision(t{3}));
        case 'POINTS'
            skipBlock(fid,3*toCount(t{2},line),vtkPrecision(t{3}));
        case 'CELLS'
            skipBlock(fid,toCount(t{3},line),'int32');
        case 'CELL_TYPES'
            skipBlock(fid,toCount(t{2},line),'int32');
        case 'CELL_DATA'
            cellCount=toCount(t{2},line); vars=readAttributes(fid,cellCount,requested); break;
        case 'POINT_DATA'
            skipAttributes(fid,toCount(t{2},line));
        otherwise
            % Other ASCII metadata is ignored until CELL_DATA is reached.
    end
end
if isempty(cellCount), error('CELL_DATA was not found in %s',filePath); end
output=struct('cellCount',cellCount,'variables',vars);
end

function vars=readAttributes(fid,n,requested)
vars=struct();
while ~feof(fid)
    line=fgetl(fid); if ~ischar(line), break; end
    line=strtrim(line); if isempty(line), continue; end
    t=strsplit(line); key=upper(t{1});
    switch key
        case 'SCALARS'
            if numel(t)<3, error('Malformed SCALARS line: %s',line); end
            name=t{2}; precision=vtkPrecision(t{3}); components=1;
            if numel(t)>=4, components=toCount(t{4},line); end
            lookup=nonemptyLine(fid);
            if ~startsWith(upper(strtrim(lookup)),'LOOKUP_TABLE'), error('Missing LOOKUP_TABLE after %s',name); end
            count=n*components;
            if ismember(name,requested)
                values=fread(fid,count,['*' precision]);
                if numel(values)~=count, error('Unexpected EOF reading %s',name); end
                if components==1, vars.(name)=double(values); else, vars.(name)=reshape(double(values),components,[]).'; end
                consumeEol(fid);
            else
                skipBlock(fid,count,precision);
            end
        case 'VECTORS'
            if numel(t)<3, error('Malformed VECTORS line: %s',line); end
            skipBlock(fid,3*n,vtkPrecision(t{3}));
        case 'FIELD'
            arrays=toCount(t{3},line);
            for a=1:arrays
                hdr=strsplit(strtrim(nonemptyLine(fid)));
                if numel(hdr)<4, error('Malformed FIELD array declaration.'); end
                skipBlock(fid,toCount(hdr{2},line)*toCount(hdr{3},line),vtkPrecision(hdr{4}));
            end
        otherwise
            % DARSim's writers use SCALARS and VECTORS for CELL_DATA.
    end
end
end

function skipAttributes(fid,n)
while ~feof(fid)
    pos=ftell(fid); line=fgetl(fid); if ~ischar(line), return; end
    line=strtrim(line); if isempty(line), continue; end
    t=strsplit(line); key=upper(t{1});
    if strcmp(key,'CELL_DATA'), fseek(fid,pos,'bof'); return; end
    switch key
        case 'SCALARS'
            p=vtkPrecision(t{3}); c=1; if numel(t)>=4, c=toCount(t{4},line); end
            nonemptyLine(fid); skipBlock(fid,n*c,p);
        case 'VECTORS'
            skipBlock(fid,3*n,vtkPrecision(t{3}));
        otherwise
    end
end
end

function skipBlock(fid,count,precision)
bytes=bytesPer(precision);
if fseek(fid,count*bytes,'cof')~=0, error('Unable to skip binary VTK block.'); end
consumeEol(fid);
end

function consumeEol(fid)
pos=ftell(fid); b=fread(fid,1,'*uint8');
if isempty(b), return; end
if b==10, return; end
if b==13
    b2=fread(fid,1,'*uint8'); if ~isempty(b2)&&b2~=10, fseek(fid,-1,'cof'); end
else
    fseek(fid,pos,'bof');
end
end

function line=nonemptyLine(fid)
while true
    line=fgetl(fid); if ~ischar(line), error('Unexpected end of VTK file.'); end
    if ~isempty(strtrim(line)), return; end
end
end

function n=toCount(text,line)
n=str2double(text);
if isnan(n)||n<0||n~=floor(n), error('Invalid count in: %s',line); end
end

function p=vtkPrecision(type)
switch lower(type)
    case {'float','float32'}, p='single';
    case {'double','float64'}, p='double';
    case {'int','int32'}, p='int32';
    case {'unsigned_int','uint32'}, p='uint32';
    case {'short','int16'}, p='int16';
    case {'unsigned_short','uint16'}, p='uint16';
    case {'char','int8'}, p='int8';
    case {'unsigned_char','uint8'}, p='uint8';
    otherwise, error('Unsupported VTK numeric type: %s',type);
end
end

function b=bytesPer(p)
switch p
    case {'double','int64','uint64'}, b=8;
    case {'single','int32','uint32'}, b=4;
    case {'int16','uint16'}, b=2;
    case {'int8','uint8'}, b=1;
    otherwise, error('Unsupported precision: %s',p);
end
end
