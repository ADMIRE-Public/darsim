%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [total_new_area] = find_all_projected_area(...
                    face_projected_ind, int_area, all_normals, eq)

loc_normals = all_normals(face_projected_ind,:);
plane_normals = repmat([eq.A eq.B eq.C],size(loc_normals,1),1);
ang_cos = dot(loc_normals,plane_normals,2);
total_new_area = abs(int_area.*ang_cos);
end
