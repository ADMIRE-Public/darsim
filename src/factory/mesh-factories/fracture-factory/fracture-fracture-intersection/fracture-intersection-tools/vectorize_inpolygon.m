%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ref] = vectorize_inpolygon(points, vol_query,...
                       poly_coords,poly_face_adj, poly_node_adj, faces_adj)
%vectorize_inpolygon checks if points are located in or out of a set of
%poluygons.
% Input:
%   points        -> nx2 matrix containing the query points
%   vol_querry    -> nx1 vector containing the indices of the polygons to
%   which a point needs to be tested
%   poly_coords   -> the coordinates of the polygons
%   poly_face_adj -> the face adjacencies of the tested polygon.
%   poly_node_adj -> the node adjacencies of the tested polygon.
%   faces_adj     -> the adjacencies of the faces that comprise the polygon
% Output:
%   ref - 1 if point in points is inside volume in vol_query
%% for test purpose only
%   poly_face_neigh -> nx2 face_neighbors, boundary elements repeat the
%                       same neighbours
%% intializing rotation matrix
% k = 1000;
tol = eps ;
teta = pi/2;
R = [cos(teta) -sin(teta); sin(teta) cos(teta)];
% defining centers 
%% initalizing face parameters
%face center -> fcenter
fcenter = 0.5*(poly_coords(faces_adj(:,1),1:2) + ...
                                        poly_coords(faces_adj(:,2),1:2));
%face_vec -> unitary vector from p2-p1 from a face p1 p2
face_vec = poly_coords(faces_adj(:,2),1:2) - ...
                                           poly_coords(faces_adj(:,1),1:2);
face_vec = face_vec./vecnorm(face_vec,2,2);
% defining normals by rotating face_vecs
normal = (R*face_vec')';
%% calculating element center skipping 0 indices
% cent_elem -> element or cell center
num_node = zeros(size(poly_node_adj,1),1);
cent_elem = zeros(size(poly_node_adj,1),2);
for index = 1:size(poly_node_adj,2)
    ref = poly_node_adj(:,index) ~=0;
    num_node(ref) = num_node(ref) + 1;
    cent_elem(ref,:) = cent_elem(ref,:) + ...
                            poly_coords(poly_node_adj(ref,index),1:2);
end
cent_elem = cent_elem ./ num_node;
% %% defining normals orientation orientation
% % vec_orient is a vector that contains that links the second neighbour to 
% % the first in poly_face_neighbour. if the a given face is a boundary
% % element the connects the negihbour centroid to the face centroid.
% % if the dot product between vec_orient and the normals is negative the
% % normals are multiplied by -1.
% vec_orient = ...
%     cent_elem(poly_face_neigh(:,2),:) - cent_elem(poly_face_neigh(:,1) ,:);
% ref = all(vec_orient == 0,2);
% vec_orient(ref,:) =  fcenter(ref,:);
% vec_orient = vec_orient./ vecnorm(vec_orient,2,2);
% ref_negative_orientation = ~(dot(normal,vec_orient,2) > 0);
% normal(ref_negative_orientation,:) = -normal(ref_negative_orientation,:);

% %% checks if point in points belong to respective vol_query
% % points_normal_list contains lists the normals of the vol_query 
% % cent_query is vector oriented the center of vol_query to the analyzed
% % points.
% % sign_matrix contain the results of the test that checks if the normals 
% % and cent_querry are aligned. 
% points_normal_list = poly_face_adj(vol_query,:);
% cent_query = points - cent_elem(vol_query,:);
% sign_matrix = zeros(size(points_normal_list));
% for index = 1:size(points_normal_list,2)
%     lref = points_normal_list(:,index) ~= 0;
%     vec_check = points(lref,:) - fcenter(points_normal_list(lref,index),:);
%     vec_check = vec_check./ vecnorm(vec_check,2,2);    
%     vec_orient = points(lref,:) - cent_query(lref,:);
%     dot_norm_point_vec = ...
%         ~(dot(normal(points_normal_list(lref,index),:), vec_check,2) > 0);
%     dot_norm_cent_query = ...
%         ~(dot(normal(points_normal_list(lref,index),:), vec_orient,2) > 0);
%     lref_aligned = dot_norm_point_vec == dot_norm_cent_query;
%     lref_pos = lref;
%     lref_pos(lref_pos) = lref_aligned;
%     lref_neg = lref;
%     lref_neg(lref_neg) = ~lref_aligned;
%     sign_matrix(lref_pos,index) = 1;
%     sign_matrix(lref_neg,index) = -1;
% end
% sign_matrix(sign_matrix == 0) = 1;
% ref = all(sign_matrix == 1,2);
%% checks if point in points belong to respective vol_query
% points_normal_list contains lists the normals of the vol_query 
% cent_query is vector oriented the center of vol_query to the analyzed
% points.
% sign_matrix contain the results of the test that checks if the normals 
% and cent_querry are aligned. 
points_normal_list = poly_face_adj(vol_query,:);
% cent_query = points - cent_elem(vol_query,:);
sign_matrix = zeros(size(points_normal_list));
for index = 1:size(points_normal_list,2)
    lref = points_normal_list(:,index) ~= 0;
    vec_check = points(lref,:) - fcenter(points_normal_list(lref,index),:);
    vec_check = vec_check./ halfvecnorml2(vec_check);    
    cent_query = cent_elem(vol_query(lref),:) - ...
                            fcenter(points_normal_list(lref,index),:);    
    cent_query = cent_query./ halfvecnorml2(cent_query);   


    %vec_orient = points(lref,:) - cent_query(lref,:);
    % dot_norm_point_vec = ...
    %     ~(dot(normal(points_normal_list(lref,index),:), vec_check,2) > 0);
    % dot_norm_cent_query = ...
    %     ~(dot(normal(points_normal_list(lref,index),:), cent_query,2) > 0);
    % 
    dot_norm_pointvec_aux = ...
            dot(normal(points_normal_list(lref,index),:), vec_check,2);
    dot_norm_point_vec_aux = ...
            dot(normal(points_normal_list(lref,index),:), cent_query,2);

    dot_norm_pointvec_aux(abs(dot_norm_pointvec_aux) < tol) = 0;
    dot_norm_point_vec_aux(abs(dot_norm_point_vec_aux) < tol) = 0;

    dot_norm_point_vec = ~(dot_norm_pointvec_aux >= 0);      
    dot_norm_cent_query = ~(dot_norm_point_vec_aux >= 0);

    lref_aligned = dot_norm_point_vec == dot_norm_cent_query;
    lref_pos = lref;
    lref_pos(lref_pos) = lref_aligned;
    lref_neg = lref;
    lref_neg(lref_neg) = ~lref_aligned;
    sign_matrix(lref_pos,index) = 1;
    sign_matrix(lref_neg,index) = -1;
end
sign_matrix(sign_matrix == 0) = 1;
ref = all(sign_matrix == 1,2);
end
