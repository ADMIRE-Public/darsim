%% Startup Routine intialize MeshCore and initialization routines
startup;
%% using gmsh parser
file_mesh = '/mnt/sda1/projects/meshes/2d/netest.msh';
A = py.gmsh_extract.GmshCore(file_mesh);

%% Using ECLIPSE File Forrmat to read Johanson Case.
%grdecl_path = "/mnt/sda1/Artur/project/Artur_DARSim/DARSim/Input/Examples/NPD5/NPD5_Z_6-10.grdecl";
grdecl_path = "/mnt/sda1/Artur/project/Artur_DARSim/DARSim/Input/Examples/Geologically_Relevant_Fields/Johansen/NPD5_Section/NPD5.grdecl";
% grdecl_file = readGRDECL(grdecl_path);
% usys   = getUnitSystem('METRIC');
% grdecl_file = convertInputUnits(grdecl_file, usys);
% G = processGRDECL(grdecl_file, 'Verbose', true);
% rock = grdecl2Rock(grdecl_file, G.cells.indexMap);
% core = MeshCore(G,"props", rock,"storeG", true);
exx = struct("permxx",1,"permyy",2);
G = grdecl_path;
core = MeshCore(G,"storeG", true, "extraProps",exx);
core.plot
%% mesh processing
% core = MeshCore(G);
