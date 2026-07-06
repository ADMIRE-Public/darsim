function [ncoords] = proj_coplanar_coords_on_xy(base_point,...
    pnormal,coords, flag)
center = base_point;
centerm = repmat(center,size(coords,1),1);
rotM = rotate_to_z_axis(pnormal);
projvec = coords - centerm;
if flag
    ncoords = (rotM*projvec')' + centerm;
else
    ncoords = ((rotM^-1)*projvec')' + centerm;
end
end