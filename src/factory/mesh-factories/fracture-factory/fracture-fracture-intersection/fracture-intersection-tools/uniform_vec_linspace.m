%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [points_2d, all_ind] = uniform_vec_linspace(xyz_start,xyz_end,n)
ind_count = [1:size(xyz_start,1)]';
dtotal = (xyz_end - xyz_start)./(n-1);
nslices = [0:(n-1)]';
ii = repelem(ind_count,size(nslices,1),1);
jj = repmat(nslices,size(dtotal,1),1);
%kk = repelem(nslices,size(dtotal,1),1);
points = xyz_start(ii,:) + jj.*[dtotal(ii,1) zeros(size(ii))]; 
all_ind = repelem(ind_count,n*n,1);
points_ind = [1:size(points,1)]';
ii = repelem(points_ind,n,1);
jj = repmat(jj,n,1);
points_2d =  points(ii,:) + jj.* [zeros(size(jj)) dtotal(all_ind,2)]; 
end
