
file_path = append(pwd,"/res_mesh.m");
[G,n] = gmshToMRST(file_path);
meshcore = MeshCore(G);
meshcore = translateMeshCore(meshcore,1);    
%core.plot
