%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [int_rock_frac_uniq, int_area] = ...
    calculate_int_area_edfm(is_in, frac_el, frac_area)
% note: one node in the fracelemt can belong simultaniously to
% two rockmesh elements (node on a shared face). Each node is therefore
% weighted by 1/(number of rock elements containing it). The area of a
% fracture element is normalized by all nodes generated in it, so nodes
% outside the rock domain carry their share of the area out of the domain.
node_mult = accumarray(double(is_in(:,2)), 1, [numel(frac_el) 1]);
node_weight = 1 ./ node_mult(is_in(:,2));
frac_point_count = accumarray(double(frac_el(:)), 1, [numel(frac_area) 1]);
int_rockfrac = [is_in(:,1)  frac_el(is_in(:,2))];
[int_rock_frac_uniq,~,ic] = unique(int_rockfrac,"rows");
h = accumarray(ic,node_weight);
int_area = h./ frac_point_count(int_rock_frac_uniq(:,2)) .* ...
    frac_area(int_rock_frac_uniq(:,2));
end
