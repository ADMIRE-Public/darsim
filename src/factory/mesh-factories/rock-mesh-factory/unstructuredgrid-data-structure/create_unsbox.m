function meshcore = create_unsbox(celldim,physdim)
G = cartGrid(celldim, physdim);
core = MeshCore(G);
meshcore = translateMeshCore(core,1);  
end