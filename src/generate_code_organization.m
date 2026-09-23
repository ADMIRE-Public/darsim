function generate_code_maps()
% Generate DARSim HTML code maps from the current src tree.
% Put this file in src and run: generate_code_maps

SRC = fileparts(mfilename('fullpath'));
maps = {
 'code-organization.html', SRC;
 fullfile('darsim-interface','darsim-interface.html'), fullfile(SRC,'darsim-interface');
 fullfile('darsim_legacy','darsim_legacy.html'), fullfile(SRC,'darsim_legacy');
 fullfile('factory','factory.html'), fullfile(SRC,'factory');
 fullfile('reader','reader.html'), fullfile(SRC,'reader');
 fullfile('simulation','simulation.html'), fullfile(SRC,'simulation');
 fullfile('tools','tools.html'), fullfile(SRC,'tools')};

mapPages = cellfun(@(x) canon(fullfile(SRC,x)), maps(:,1), 'UniformOutput',false);

for k=1:size(maps,1)
    page=fullfile(SRC,maps{k,1}); root=maps{k,2};
    if ~isfolder(fileparts(page)), mkdir(fileparts(page)); end
    intro=extractIntro(page);
    if same(root,SRC), heading='DARSim Organization'; title='DARSim Code Organization';
    else, [~,heading]=fileparts(root); title=[heading ' — DARSim Code Map']; end
    body=treeHTML(root,root,fileparts(page),SRC,maps,mapPages);
    html=sprintf(['<!DOCTYPE html>\n<html lang="en"><head><meta charset="UTF-8">' ...
      '<meta name="viewport" content="width=device-width, initial-scale=1">' ...
      '<meta name="generator" content="generate_code_maps.m"><title>%s</title>' ...
      '<style>%s</style></head><body><main class="container">' ...
      '<nav class="breadcrumbs">%s</nav><h1>%s</h1>' ...
      '<section class="comment-box">%s</section>' ...
      '<div class="toolbar"><input id="filter" type="search" placeholder="Filter files or paths…">' ...
      '<button id="expand">Expand all</button><button id="collapse">Collapse all</button>' ...
      '<span id="status" class="status"></span></div>' ...
      '<div id="code-tree" class="tree">%s</div>' ...
      '<footer>Generated from the current <code>src</code> directory by <code>generate_code_maps.m</code>.</footer>' ...
      '</main><script>%s</script></body></html>\n'], ...
      esc(title),CSS(),crumbs(page,root,SRC),esc(heading),intro,body,JS());
    fid=fopen(page,'w','n','UTF-8'); assert(fid>=0,'Cannot write %s',page);
    fwrite(fid,html,'char'); fclose(fid);
end

broken={};
for k=1:size(maps,1)
    page=fullfile(SRC,maps{k,1}); b=validate(page);
    if ~isempty(b), broken(end+1,:)={maps{k,1},b}; end %#ok<AGROW>
end
fprintf('Generated %d code maps.\n',size(maps,1));
if ~isempty(broken)
    for k=1:size(broken,1), fprintf('BROKEN %s: %s\n',broken{k,1},strjoin(broken{k,2},', ')); end
    error('One or more generated links are broken.');
end
fprintf('Validated all generated file and navigation links.\n');
end

function s=treeHTML(d,root,pageDir,SRC,maps,mapPages)
E=entries(d,mapPages); out={'<ul>'};
for i=1:numel(E)
 p=fullfile(d,E(i).name); name=esc(E(i).name);
 if E(i).isdir
   if same(root,SRC)
     dest='';
     for j=2:size(maps,1)
       [pd,~,~]=fileparts(maps{j,1}); [~,dn]=fileparts(pd);
       if strcmp(dn,E(i).name), dest=fullfile(SRC,maps{j,1}); break; end
     end
     if ~isempty(dest) && isfile(dest)
       out{end+1}=sprintf('<li class="dir-item" data-path="%s"><a class="folder" href="%s">📁 %s</a></li>',name,href(dest,pageDir),name); %#ok<AGROW>
     else
       out{end+1}=sprintf('<li class="dir-item" data-path="%s"><span class="folder">📁 %s</span></li>',name,name); %#ok<AGROW>
     end
   else
     rel=esc(strrep(relpath(p,root),filesep,'/'));
     if ~isempty(entries(p,mapPages))
       out{end+1}=sprintf('<li class="dir-item" data-path="%s"><details><summary><span class="folder">📁 %s</span></summary>%s</details></li>',rel,name,treeHTML(p,root,pageDir,SRC,maps,mapPages)); %#ok<AGROW>
     else
       out{end+1}=sprintf('<li class="dir-item" data-path="%s"><span class="folder">📁 %s</span></li>',rel,name); %#ok<AGROW>
     end
   end
 else
   rel=esc(strrep(relpath(p,root),filesep,'/'));
   out{end+1}=sprintf('<li class="file-item" data-path="%s"><a class="file" href="%s">%s</a></li>',rel,href(p,pageDir),name); %#ok<AGROW>
 end
end
out{end+1}='</ul>'; s=strjoin(out,'');
end

function E=entries(folder,mapPages)
E=dir(folder); E=E(~ismember({E.name},{'.','..'})); keep=true(size(E));
for i=1:numel(E)
 n=E(i).name; [~,~,x]=fileparts(n);
 bad=ismember(n,{'.git','.idea','.vscode','__pycache__','.pytest_cache','.mypy_cache','.DS_Store','Thumbs.db','generate_code_maps.m','generate_code_maps.py'}) || ...
     ismember(lower(x),{'.pyc','.pyo','.swp','.tmp','.bak'}) || any(strcmp(canon(fullfile(folder,n)),mapPages));
 keep(i)=~bad;
end
E=E(keep); 
if isempty(E), return; end
T=table(~[E.isdir]',lower(string({E.name}))','VariableNames',{'file','name'});
[~,ix]=sortrows(T,{'file','name'}); E=E(ix);
end

function s=extractIntro(page)
if isfile(page)
 t=fileread(page); q=regexp(t,'<div\s+class=["'']comment-box["'']\s*>(.*?)</div>','tokens','once');
 if ~isempty(q), s=strtrim(q{1}); return; end
end
[~,n,e]=fileparts(page);
if strcmp([n e],'code-organization.html')
 s=['<p>Here, you can explore the structure of the files and folders within this repository. The root folder is the <code>src</code> folder. Click a file to inspect it.</p><br>' ...
 '<p>The primary routines in this folder are <code>DARSim.m</code> and <code>startup_darsim.m</code>. <code>DARSim.m</code> initiates the simulation process, while <code>startup_darsim.m</code> contains routines that assist with initialization.</p>'];
else, s='<p>Browse the current source-code structure below.</p>'; end
end

function s=crumbs(page,root,SRC)
if same(root,SRC), s='<span>src</span>';
else, [~,n]=fileparts(root); s=sprintf('<a href="%s">src</a> / <span>%s</span>',href(fullfile(SRC,'code-organization.html'),fileparts(page)),esc(n)); end
end

function h=href(target,base)
p=strsplit(strrep(relpath(target,base),filesep,'/'),'/');
for i=1:numel(p), p{i}=urlpart(p{i}); end
h=strjoin(p,'/');
end

function s=urlpart(s)
b=unicode2native(s,'UTF-8'); safe='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-._~'; z='';
for x=b
 if x<128 && contains(safe,char(x)), z=[z char(x)]; else, z=[z sprintf('%%%02X',x)]; end %#ok<AGROW>
end
s=z;
end

function b=validate(page)
t=fileread(page); q=regexp(t,'href="([^"]+)"','tokens'); b={};
for i=1:numel(q)
 h=q{i}{1}; if startsWith(h,{'http:','https:','#','mailto:'}), continue; end
 h=urldecode(h); p=fullfile(fileparts(page),strrep(h,'/',filesep));
 if ~(isfile(p)||isfolder(p)), b{end+1}=q{i}{1}; end %#ok<AGROW>
end
end

function s=urldecode(s)
q=regexp(s,'%([0-9A-Fa-f]{2})','tokens');
while ~isempty(q), x=q{1}{1}; s=regexprep(s,['%' x],char(hex2dec(x)),'once'); q=regexp(s,'%([0-9A-Fa-f]{2})','tokens'); end
end

function r=relpath(target,base)
t=strsplit(strrep(canon(target),filesep,'/'),'/'); b=strsplit(strrep(canon(base),filesep,'/'),'/'); k=0;
for i=1:min(numel(t),numel(b)), if peq(t{i},b{i}), k=i; else, break; end, end
p=[repmat({'..'},1,numel(b)-k),t(k+1:end)]; if isempty(p), r='.'; else, r=strjoin(p,filesep); end
end
function tf=peq(a,b), if ispc, tf=strcmpi(a,b); else, tf=strcmp(a,b); end, end
function tf=same(a,b), tf=peq(canon(a),canon(b)); end
function p=canon(p), try, p=char(java.io.File(p).getCanonicalPath()); catch, p=char(string(p)); end, end
function s=esc(s), s=strrep(char(s),'&','&amp;'); s=strrep(s,'<','&lt;'); s=strrep(s,'>','&gt;'); s=strrep(s,'"','&quot;'); s=strrep(s,'''','&#39;'); end

function s=CSS()
s=[':root{color-scheme:light dark;--bg:#fff;--fg:#202124;--muted:#5f6368;--panel:#f5f7f8;--line:#d8dee4;--link:#075cc7;--hover:#fff3a3}' ...
'@media(prefers-color-scheme:dark){:root{--bg:#16181b;--fg:#e8eaed;--muted:#aeb4bb;--panel:#22262b;--line:#3a4148;--link:#8ab4f8;--hover:#554d16}}' ...
'*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}.container{max-width:1400px;margin:auto;padding:24px}' ...
'a{color:var(--link);text-decoration:none}a:hover{text-decoration:underline;background:var(--hover)}.breadcrumbs{margin-bottom:14px;color:var(--muted)}h1{font-family:system-ui,sans-serif}' ...
'.comment-box{padding:18px 20px;background:var(--panel);border:1px solid var(--line);border-radius:8px;line-height:1.45;margin-bottom:18px}.toolbar{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin:14px 0}' ...
'#filter{flex:1 1 320px;max-width:620px;padding:9px 11px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--fg);font:inherit}button{padding:8px 10px;border:1px solid var(--line);border-radius:6px;background:var(--panel);color:var(--fg);cursor:pointer}' ...
'.tree,.tree ul{list-style:none;margin:0;padding-left:0}.tree ul{padding-left:22px;border-left:1px dotted var(--line);margin-left:7px}.tree li{line-height:1.65;white-space:nowrap}details>summary{cursor:pointer;list-style:none}details>summary::before{content:"▸ ";color:var(--muted)}details[open]>summary::before{content:"▾ "}.folder{font-weight:650}.file::before{content:"· ";color:var(--muted)}.match-hidden{display:none!important}.status,footer{color:var(--muted);font-family:system-ui,sans-serif}footer{margin-top:28px;padding-top:14px;border-top:1px solid var(--line);font-size:.85rem}'];
end

function s=JS()
s=['const filter=document.getElementById("filter"),tree=document.getElementById("code-tree"),status=document.getElementById("status");' ...
'function applyFilter(){const q=filter.value.trim().toLowerCase();let n=0;const a=[...tree.querySelectorAll("li")];if(!q){a.forEach(x=>x.classList.remove("match-hidden"));tree.querySelectorAll("details").forEach(x=>x.open=false);status.textContent="";return}a.reverse().forEach(x=>{const own=(x.dataset.path||"").toLowerCase().includes(q),child=!!x.querySelector("li:not(.match-hidden)"),show=own||child;x.classList.toggle("match-hidden",!show);const d=x.querySelector(":scope > details");if(d&&show)d.open=true;if(show&&x.classList.contains("file-item"))n++});status.textContent=n+(n===1?" matching file":" matching files")}' ...
'filter.addEventListener("input",applyFilter);document.getElementById("expand").addEventListener("click",()=>tree.querySelectorAll("details").forEach(x=>x.open=true));document.getElementById("collapse").addEventListener("click",()=>tree.querySelectorAll("details").forEach(x=>x.open=false));'];
end
