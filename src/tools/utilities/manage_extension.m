function [meshcore] = manage_extension(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[fpath,fname,ext] = fileparts(filename);
switch ext
    case ".m"
        run(filename)
        if ~(exist("meshcore") == 1)
            error("Error: %s script needs to define a meshcore variable", ...
                strjoin(fname, ext))
        end
    case ".msh"
        fname = char(fname);
        filename_gmsh = fullfile(fpath, [fname, '.m']);
        if isfile(filename_gmsh)
            if isfile(fullfile(fpath, [fname, '.mat']))
                load(fullfile(fpath, [fname, '.mat']),'Gmshcore')
            else
                Gmshcore = gmshToMRST(filename_gmsh);
                save(fullfile(fpath, [fname, '.mat']),"Gmshcore")
            end
            core = MeshCore(Gmshcore);
            meshcore = translateMeshCore(core,1);    
        else
            disp("Gmsh file in matlab format not found. Converting msh file.")
            disp("Python and MATLAB integration still unstable. MATLAB may crash.")
            pe = pyenv;
            if ~isempty(pe.Version)
                %clear classes
                command_string = strjoin(["msh-converter.py ", filename]);
                pyrunfile(command_string)
                terminate(pe)
            else
                error("Error: MATLAB cannot access Python.")
            end
        end

        
    otherwise
        error("Error: File Format not yet available")
end
end