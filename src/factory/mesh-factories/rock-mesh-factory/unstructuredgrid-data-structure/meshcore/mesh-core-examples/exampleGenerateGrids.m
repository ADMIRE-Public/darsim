%% Set of examples that explain how MRST can be used to generate the meshes
%

%% Startup Routine intialize MeshCore andinitialization routine
startup;

%% Cartesian Grid Intialization
G = cartGrid([30,20,4]);
G = removeCells(G, [2:40]);
core = MeshCore(G);
core.plot

%% Simple Cornerpoint Grid Initialization
grdecl = simpleGrdecl([20, 20, 5], ...
@(x) .05*(sin(2*pi*x)-1.5));
G = processGRDECL(grdecl, 'Verbose', true);
core = MeshCore(G);
core.plot
%% Extrude Triangular Initalization
% Triangular grids
% Generate a 2D triangle grid from the same data using the xy-coordinates
pts = rand(20,3).*repmat([10,10,1], 20, 1);
pts_2d = pts(:,1:2);
% Triangulate
T = delaunay(pts_2d(:,1), pts_2d(:,2));
G_tri = triangleGrid(pts_2d, T);
% Extruded grid
% If we have an interesting 2D grid, we can easily extend it to 3D by using
% makeLayeredGrid. We will extend the previous triangle grid to 3 layers
G = makeLayeredGrid(G_tri, 3);
core = MeshCore(G);
core.plot
%% Voronoi Grid Intialization
% Triangular grid in 2D
[x,y] = meshgrid(0:10,0:8);
x(2:8,2:10) = x(2:8,2:10) + .25*randn(7,9);
y(2:8,2:10) = y(2:8,2:10) + .25*randn(7,9);
G = triangleGrid([x(:) y(:)]);
h=plotGrid(G,'FaceColor','none');
% Voronoi diagram
V = pebi(G);
hold on, plot(x(:),y(:),'k.','MarkerSize',12); hold off
plotGrid(V,'FaceColor','none','LineWidth',2);
delete(h)
G = makeLayeredGrid(V, 3);
core = MeshCore(G);
core.plot
%% Using ECLIPSE File Forrmat to read Johanson Case.
%grdecl_path = "/mnt/sda1/Artur/project/Artur_DARSim/DARSim/Input/Examples/NPD5/NPD5_Z_6-10.grdecl";
grdecl_path = "/mnt/sda1/Artur/project/Artur_DARSim/DARSim/Input/Examples/Geologically_Relevant_Fields/Johansen/NPD5_Section/NPD5.grdecl";
grdecl_file = readGRDECL(grdecl_path);
usys   = getUnitSystem('METRIC');
grdecl_file = convertInputUnits(grdecl_file, usys);
G = processGRDECL(grdecl_file, 'Verbose', true);
rock = grdecl2Rock(grdecl_file, G.cells.indexMap);
core = MeshCore(G);
core.plot
%% mesh processing
core = MeshCore(G);
