%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function is_in = in_fracture_bounding_box(plane_geo,fp1, fp8)
[~,fp2,fp3,fp4,fp5,fp6,fp7,~] = generate_bound_box_vertices(fp1,fp8);

%all_points = [fp1; fp2; fp3; fp4; fp5; fp6; fp7; fp8];
l_edges = cell(12,1);
l_edges{1} = [fp1 fp2];
l_edges{2} = [fp2 fp3];
l_edges{3} = [fp3 fp4];
l_edges{4} = [fp4 fp1];

l_edges{5} = [fp5 fp6];
l_edges{6} = [fp6 fp7];
l_edges{7} = [fp7 fp8];
l_edges{8} = [fp8 fp1];

l_edges{9} = [fp1 fp5];
l_edges{10} = [fp2 fp6];
l_edges{11}= [fp3 fp7];
l_edges{12} = [fp4 fp8];


% tag = logical([0 0 0; ...
%                    0 0 1; ...
%                    0 1 0; ...
%                    0 1 1; ...
%                    1 0 0; ...
%                    1 0 1; ...
%                    1 1 0;...
%                    1 1 1]);
is_in =  ...
  in_planar_bounding_box_updated(plane_geo.coords, plane_geo.normal, l_edges);
% 
% frac_1 = min(plane_geo.coords);
% frac_8 = max(plane_geo.coords);
% [~,frac_2,frac_3,frac_4,frac_5,frac_6,frac_7,~] = generate_bound_box_vertices(frac_1,frac_8);
% all_points = [frac_1; frac_2; frac_3; frac_4; frac_5; frac_6; frac_7; frac_8];
% 
% ref_t = false(size(fp1,1),1);
% for index = 1:size(fp1,1)
%     ref = in_between(all_points, fp1(index,:), fp8(index,:), 10^-8);
%     1
% end



%boundbox = any(reshape(all(is_in,2),[],8),2);
end
