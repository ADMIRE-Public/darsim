%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avg_d_frac = create_frac_mesh_avg_dist(frac_mesh_int, ...
                                               fracmesh, geo_intersection)
avg_d_frac = cell(size(frac_mesh_int));
for ii = 1:size(frac_mesh_int,1)
        for jj = ii+1:size(frac_mesh_int,2)
            if ~isempty(frac_mesh_int{ii,jj})
                fract_intersection = cell2mat(...
                   geo_intersection.get_boundary_line_intersection(ii,jj));                   
                auxmat = frac_mesh_int{ii,jj};
                p1 = fracmesh(ii).elements.centroid(auxmat(:,1),:);
                p2 = fracmesh(jj).elements.centroid(auxmat(:,2),:);
                [d1,d2] = dist_point_to_line(p1, p2, ...
                        fract_intersection);      
                avg_d_frac{ii,jj} = [d1,d2];
            end
        end
end

end

function [d1,d2] = dist_point_to_line(p1, p2, line)
    %% https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
%%
x0 = [p1;p2];
x1 = repmat(line(1,:),size(x0,1),1);
x2 = repmat(line(2,:),size(x0,1),1);
A = vecnorm(cross((x0-x1), (x0-x2)),2,2);
B = norm(p1-p2,2);
d = A/B;
d1 = d(1:size(p1,1));
d2 = d(size(p1,1)+1:end);
end
