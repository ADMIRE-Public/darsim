%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [int_rock_frac_uniq, int_area] = ...
    calculate_int_area_edfm(is_in, frac_el, frac_area)
% note: one node in the fracelemt can belong simultaniously to
% two rockmesh element. this why the percentage of the area distributed
% between these elements can be higher than 100%. to solve this instead
% of counting the original frac_el element, we will count the node
% distribution on the frac rock matrix intersection.
%frac_point_count = accumarray(frac_el,1);
frac_point_count = accumarray(frac_el(is_in(:,2)),1);
int_rockfrac = [is_in(:,1)  frac_el(is_in(:,2))];
[int_rockfrac_sorted] = sortrows(int_rockfrac);
[int_rock_frac_uniq,~,ic] = unique(int_rockfrac_sorted,"rows");
h = accumarray(ic,1);
int_area = h./ frac_point_count(int_rock_frac_uniq(:,2)) .* ...
    frac_area(int_rock_frac_uniq(:,2));
end
