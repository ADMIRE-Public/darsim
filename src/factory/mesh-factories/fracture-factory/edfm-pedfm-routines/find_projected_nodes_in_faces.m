%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out_ref] = find_projected_nodes_in_faces(points, normals, ...
                            faces_center, face_to_points, rockmesh)
projplane = zeros(size(faces_center,1),3);
projplane(:,3) = 1;
R1d =  rotate_to_z_axis_vec(normals, projplane);
projvec = points - faces_center;
coords_2d = matvecprod1dvec(R1d,projvec) + faces_center;
out_ref = false(size(face_to_points));
[faces_edges] = create_matedgezero(face_to_points);
for index = 1:(size(faces_edges,2)-1)  
    %index
    %warning("pEDFM debug 1")
    edge = faces_edges(:,index:index+1);
    ref = ~any(edge == 0,2);
    if any(ref)
        c1 = rockmesh.vertex.coord(edge(ref,1),:);
        projvec = c1 - faces_center(ref,:);
        c1_2d = matvecprod1dvec(R1d,projvec) + faces_center(ref,:);
        c2 = rockmesh.vertex.coord(edge(ref,2),:);
        projvec = c2 - faces_center(ref,:);
        c2_2d = matvecprod1dvec(R1d,projvec) + faces_center(ref,:);
        lref = is_above_2d(c1_2d,c2_2d,coords_2d,faces_center);
        out_ref(ref,index) = lref;
    end
end
out_ref(face_to_points == 0) = true;
out_ref = all(out_ref,2);   
end


