%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frac_int_avg_d] = calculate_avg_frac_d(frac_int_avg_d, ...
    test_points, frac_el, frac_avg_d_ind, frac_int, plane_int,index)
matref = ~cellfun(@isempty,frac_int);
[ii,jj] = find(matref);
ref1 = ii == index;
ref2 = jj == index;
ref = ref1|ref2;
ii = ii(ref);
jj = jj(ref);
for counter = 1:size(ii,1)
   local_ii = ii(counter);
   local_jj = jj(counter);
   ind_vol = frac_avg_d_ind{local_ii,local_jj}{~(local_ii == index)+1};
   ref = ismember(frac_el,ind_vol);
   local_frac_el = frac_el(ref);
   fracline = ...
       cell2mat(plane_int.get_boundary_line_intersection(local_ii,local_jj));
   lavg_dist = dist_point_to_line(test_points(ref,:), fracline);
   sumdist = accumarray(local_frac_el,lavg_dist);
   ndist = accumarray(local_frac_el,1);

   avg_dist = sumdist(ind_vol)./(ndist(ind_vol));
   if isempty(frac_int_avg_d{local_ii,local_jj})
       frac_int_avg_d{local_ii,local_jj} = cell(2,1);
   end
   auxcell = frac_int_avg_d{local_ii,local_jj};
   auxcell{~(local_ii == index)+1} = avg_dist;
   frac_int_avg_d{local_ii,local_jj} = auxcell;
end
end

function [d] = dist_point_to_line(p,  line)
    %% https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
%%
x0 = p;
x1 = repmat(line(1,:),size(x0,1),1);
x2 = repmat(line(2,:),size(x0,1),1);
A = vecnorm(cross((x0-x1), (x0-x2)),2,2);
B = norm(line(1,:) - line(2,:),2);
d = A./B;
end
