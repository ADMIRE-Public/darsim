function startup
warning('off','MATLAB:dispatcher:nameConflict')
%% adding matlab path
addpath(genpath(pwd))
%% adding python path
folder = fullfile("py-code");
if count(py.sys.path,folder) == 0
    insert(py.sys.path,int32(0),folder);
end
disp("MeshCore Startup successfully initialized.")
end