%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ind_struct] = create_avg_d_frac_structure(matcell)
ind_struct = cell(size(matcell,1));
for ii = 1:size(matcell,1)
    for jj = 1:size(matcell,2)
        if ~isempty(matcell{ii,jj})
            auxmat = matcell{ii,jj};
            col1 = unique(auxmat(:,1));
            col2 = unique(auxmat(:,2));            
            ind_struct{ii,jj} = {col1,col2};
        end
    end
end
end
