%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef translate_newDarsim_to_oldDarsim_old
    properties
        num_frac_el
        num_rock_el
        frac_ind
        fracmesh
        frac_num
        frac_global
        frac_local
        CrossConnections
    end
    methods
        function obj = translate_newDarsim_to_oldDarsim(new_struct, rockmesh, fracmesh)
            obj.fracmesh = fracmesh;
            obj.num_rock_el = size(rockmesh.element.vertices,1);
            [obj.num_frac_el, obj.frac_ind] = ...
                                        obj.define_frac_num(fracmesh);
            [obj.frac_num,obj.frac_global,obj.frac_local] = ....
                                obj.create_fracmesh_dictionary();

            CrossConnections = generate_empty_cross_connections(obj);
            [CrossConnections] = obj.fillCrossConnections(...
                                       CrossConnections, new_struct);
            obj.CrossConnections = fill_zeros(CrossConnections);
            
        end


        function [frac_num,frac_global,frac_local] = ...
                                create_fracmesh_dictionary(obj)
            frac_num = repelem([1:length(obj.frac_ind)]', obj.frac_ind);
            frac_global = [1:length(frac_num)]';
            frac_local = zeros(sum(obj.num_frac_el),1);
            acum = 0;
            init= 1 ;
            for index = obj.num_frac_el'
                frac_local((init+acum):(index+acum)) = [1:index]';
                acum = acum + index;
            end            
        end
        function CrossConnections = ....
                fillCrossConnections(obj,CrossConnections, new_struct)
            %load Discretizations.mat
            loc_sum = [0; obj.frac_ind];
            for index = 1:length(new_struct.fracrock_edfm)
                local_edfm = new_struct.fracrock_edfm{index};
                [frac_ind_sorted,ref] = sort(local_edfm.frac_ind);
                CI_sorted = local_edfm.CI(ref);
                vol_ind_sorted = local_edfm.vol_ind(ref);
                frac_ind_sorted = frac_ind_sorted + loc_sum(index);
                for local_index = 1:length(frac_ind_sorted)
                    frac_el = frac_ind_sorted(local_index);
                    CrossConnections(frac_el).Cells = [...
                        CrossConnections(frac_el).Cells; ...
                        vol_ind_sorted(local_index)];
                    CrossConnections(frac_el).CI = [...
                        CrossConnections(frac_el).CI;...
                        CI_sorted(local_index)];                   
                end
            end
            if new_struct.pedfm_on
                for index = 1: length(new_struct.fracrock_pedfm)
                    local_pedfm = new_struct.fracrock_pedfm{index};
                    [frac_ind_sorted,ref] = sort(local_pedfm.frac_ind);
                    CI_sorted = local_pedfm.CI(ref);
                    vol_ind_sorted = local_pedfm.vol_ind(ref);
                    frac_ind_sorted = frac_ind_sorted + loc_sum(index);
                    for local_index = 1:length(frac_ind_sorted)
                        frac_el = frac_ind_sorted(local_index);
                        CrossConnections(frac_el).Cells = [...
                            CrossConnections(frac_el).Cells; ...
                            vol_ind_sorted(local_index)];
                        CrossConnections(frac_el).CI = [...
                            CrossConnections(frac_el).CI;...
                            CI_sorted(local_index)];                   
                    end
                end
            end
            nfrac = size(new_struct.fracfrac,1);
            for ii = 1:nfrac
                for jj = 1:nfrac
                    if ~isempty(new_struct.fracfrac{ii,jj})
                        physical = new_struct.fracfrac{ii,jj};
                        ind1 = physical.ind1;
                        cross_ind1 = ind1 + loc_sum(ii);
                        ind2 = physical.ind2;
                        cross_ind2 = ind2 + loc_sum(jj);
                        [~,ref] = sortrows([cross_ind1 cross_ind2]);
                        cross_ind1 = cross_ind1(ref);
                        cross_ind2 = cross_ind2(ref);
                        CI = physical.CI(ref);
                        for local_index = 1:length(cross_ind1)
                            frac_el = cross_ind1(local_index);
                            frac_el2 = cross_ind2(local_index) + ...
                                                    obj.num_rock_el;
                            CrossConnections(frac_el).Cells = [...
                                CrossConnections(frac_el).Cells; frac_el2];
                            CrossConnections(frac_el).CI = [...
                                CrossConnections(frac_el).CI;...
                                CI(local_index)]; 
                        end
                        for local_index = 1:length(cross_ind1)
                            frac_el = cross_ind2(local_index);
                            frac_el2 = cross_ind1(local_index) + ...
                                                    obj.num_rock_el;
                            CrossConnections(frac_el).Cells = [...
                                CrossConnections(frac_el).Cells; frac_el2];
                            CrossConnections(frac_el).CI = [...
                                CrossConnections(frac_el).CI;...
                                CI(local_index)]; 
                        end


                    end
                end
            end
        end
        function [CrossConnections] = generate_empty_cross_connections(obj)
            lcell = cell(1);
            CrossConnections = struct("Cells",lcell, "CI", [], ...
                "T_Geo", [], "T_Geo_Cond",[], "UpWind",[], "U_Geo", []);
            CrossConnections(sum(obj.num_frac_el),1).CI = [];
        end
        
        function [frac_num_el,frac_el_ind] = define_frac_num(obj, fracmesh)
            num_frac =length(fracmesh);
            frac_num_el = zeros(num_frac,1);
            
            for index = 1:num_frac
                frac_num_el(index) = length(fracmesh(index).elements.area);                
            end
            frac_el_ind = cumsum(frac_num_el);
        end

    end
end


function [CrossConnections] = fill_zeros(CrossConnections)
for index = 1:length(CrossConnections)
    num_cross = size(CrossConnections(index).Cells,1);
    CrossConnections(index).T_Geo = zeros(num_cross,1);
    CrossConnections(index).T_Geo_Cond = ...
        zeros(num_cross,1);
    CrossConnections(index).UpWind = ...
        zeros(num_cross,2);
    CrossConnections(index).U_Geo = ...
        zeros(num_cross,2);
end
end
