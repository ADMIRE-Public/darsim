%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num_vol_int,d1,d2] =  ...
                           compute_mesh_interval(mesh,num_int,p_1,p_2, ...
                                                       fp1,fp2,vol1,vol2)
%% elements faces is organzies inner + bound
%% creating volume intersect nodes (vp1 vp2)
vp1 = sparse(size(mesh.elements.vertices,1),3);
vp2 = sparse(size(mesh.elements.vertices,1),3);
%% creating entities to classify and analyze volumes
vol_count = zeros(size(mesh.elements.faces,1),1);
pos_col = false(size(vol_count));
for index = 1:size(mesh.elements.faces,2)
    % ref = (mesh.elements.faces(:,index) ~= 0)  & ...
    %        num_int(mesh.elements.faces(:,index)) ~=0;
    ref1 = (mesh.elements.faces(:,index) ~= 0)  & ...
           num_int(mesh.elements.faces(:,index)) == 1;      
    ref2 = (mesh.elements.faces(:,index) ~= 0)  & ...
           num_int(mesh.elements.faces(:,index)) == 2;   
    ref1p1 = ref1 & ~pos_col;
    ref1p2 = ref1 & pos_col;
    pos_col(ref1p1) = true;
    vp1(ref1p1,:) = p_1(mesh.elements.faces(ref1p1,index),:);
    vp2(ref1p2,:) = p_1(mesh.elements.faces(ref1p2,index),:);      
    vp1(ref2,:) = p_1(mesh.elements.faces(ref2,index),:);
    vp2(ref2,:) = p_2(mesh.elements.faces(ref2,index),:);    
end
%% creating volume int count
num_vol_int = mesh.elements.faces;
ref = num_vol_int~=0;
num_vol_int(ref) = num_int(mesh.elements.faces(ref));
num_vol_int = sum(num_vol_int,2);
%% adding vol1-fp1 and vol2-fp2
if ~isempty(vol1) && any(num_vol_int(vol1) == 1)
    ref = num_vol_int(vol1) == 1;
    vp2(vol1(ref),:) = fp1;
    num_vol_int(vol1(ref)) = 2;
end
if ~isempty(vol2) && any(num_vol_int(vol2) == 1)
    ref = num_vol_int(vol2) == 1;
    vp2(vol2(ref),:) = fp2;
    num_vol_int(vol2(ref)) = 2;
end
%% creating distance intervals. distances are calculates from fp1
double_crossed_vol = num_vol_int == 2;
d1 = sparse(size(double_crossed_vol,1),1);
d2 = d1;
d1(double_crossed_vol) = vecnorm(vp1(double_crossed_vol,:) - ...
                repmat(fp1,sum(double_crossed_vol),1),2,2);
d2(double_crossed_vol) = vecnorm(vp2(double_crossed_vol,:) - ...
                repmat(fp1,sum(double_crossed_vol),1),2,2);

end
