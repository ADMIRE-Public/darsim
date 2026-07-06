%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FracMeshIntersection
    properties
        meshFracInterval
    end
    properties (Access = private)
        intersectWith
        meshLineIntersection
    end
    methods
        function obj = FracMeshIntersection(geo_intersection, ...
                                                           geometry,meshes)
            obj.intersectWith =  geo_intersection.intersectWith;
            frac_mesh_int_distances = ...
                   obj.calculate_volume_and_edges_fracline_intersect(...
                                    geo_intersection, geometry,meshes);
            obj.meshFracInterval = ...
                            obj.define_fracmesh_fracmesh_intersection(...
                                                frac_mesh_int_distances);

        end
    end
    methods (Access = private)
        function frac_mesh_intersection = ...
             calculate_volume_and_edges_fracline_intersect(obj, ...
                                        geo_intersection,geometry,meshes)
            frac_mesh_intersection = cell(size(obj.intersectWith));
            for index = find(obj.intersectWith)'
                [frac_mesh,intersection_line] = ...
                                ind2sub(size(obj.intersectWith),index);
                bound_intersection = cell2mat(geo_intersection. ...
                          get_boundary_line_intersection(frac_mesh, ...
                                                intersection_line));                
                fracture_p1 = bound_intersection(1,:);
                fracture_p2 = bound_intersection(2,:);

                [num_int, int_p1,int_p2,v1,v2]  = ...
                                calculate_vol_edge_intersection( ...
                                 fracture_p1, fracture_p2,...
                                 geometry(frac_mesh).normal,...
                                 meshes(frac_mesh));
                dist_frac_int = struct();
                
                [dist_frac_int.ref,dist_frac_int.d1, dist_frac_int.d2]= ...
                     compute_mesh_interval(meshes(frac_mesh), ...
                                      num_int, int_p1,int_p2, ...
                                      fracture_p1,fracture_p2,v1,v2);
                       
                frac_mesh_intersection{frac_mesh,intersection_line} = ...
                                                        dist_frac_int;
            end
        end
        function meshFracInterval = ...
              define_fracmesh_fracmesh_intersection(obj,frac_mesh_distance)
            meshFracInterval = cell(size(obj.intersectWith));
            [fmesh_ii, fmesh_jj] = find(obj.intersectWith);
            auxmat = unique(sort([fmesh_ii fmesh_jj],2),"rows");
            fmesh_ii = auxmat(:,1);
            fmesh_jj = auxmat(:,2);
            for index = 1:size(fmesh_jj,1)
                ii = fmesh_ii(index);
                jj = fmesh_jj(index);
                dist_ii = frac_mesh_distance{ii,jj};
                dist_jj = frac_mesh_distance{jj,ii};

                d_matrix =  compute_interval_intersection( ...
                    dist_ii.ref, dist_ii.d1, dist_ii.d2, ...
                    dist_jj.ref, dist_jj.d1, dist_jj.d2);
                meshFracInterval{ii,jj} = d_matrix;
            end
        end
    end
end

function d_matrix = compute_interval_intersection(f1ref,f1d1,f1d2,...
                                                        f2ref,f2d1,f2d2)
    ref = (f1ref == 2) & all((f1d1 ~= f1d2),2); 
    volI1 = find(ref);
    I1 = full([f1d1(ref) f1d2(ref)]);
    I1 = sort(I1,2);
    
    ref = f2ref == 2 & all((f2d1 ~= f2d2),2);
    volI2 = find(ref);
    I2 = full([f2d1(ref) f2d2(ref)]);
    I2 = sort(I2,2);
    refAI1 = repmat(volI1, size(volI2,1),1);
    indI1 = [1:size(volI1,1)]';
    repindI1 = repmat(indI1, size(volI2,1),1);
    
    refAI2 = repelem(volI2, size(volI1,1),1);
    indI2 = [1:size(volI2,1)]';
    repindI2 = repelem(indI2, size(volI1,1),1);
    interval_int = interval_intersection(I1(repindI1,:),I2(repindI2,:));
    d = abs(interval_int(:,2) -  interval_int(:,1));
    ref = (d~=0);
    d_matrix = [refAI1(ref),refAI2(ref), d(ref)];
end


