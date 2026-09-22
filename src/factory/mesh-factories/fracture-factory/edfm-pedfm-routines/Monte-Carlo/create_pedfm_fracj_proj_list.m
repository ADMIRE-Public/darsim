%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pedfm_vol_frac, area] = create_pedfm_fracj_proj_list( ...
    pedfm_face, pedfm_volume, original_vol, rock_frac_ind, int_area)
[rock_frac_ind,ref] = sortrows([rock_frac_ind(:,2) rock_frac_ind(:,1)]);
int_area = int_area(ref);
[~,~,cref] = unique(pedfm_volume);
num_elem = accumarray(cref,1);
int_area = int_area(ref);
end
