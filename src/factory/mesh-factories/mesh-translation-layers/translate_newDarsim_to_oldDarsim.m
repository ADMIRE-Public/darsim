%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% Translating new DARSim fracture connectivity data to the old DARSim
% cross connection format.
%
% Input:
%   new_struct:     struct
%                   Contains fracture-rock and fracture-fracture
%                   connectivity information from the new DARSim format
%   rockmesh:       struct
%                   Contains relevant rock mesh information
%   fracmesh:       struct
%                   Contains relevant fracture mesh information
%
% Output:
%   CrossConnections: struct
%                     Contains cross connection information for each
%                     fracture cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef translate_newDarsim_to_oldDarsim

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

        function obj = ...
                translate_newDarsim_to_oldDarsim(new_struct, ...
                                                 rockmesh, ...
                                                 fracmesh)

            % Store fracture mesh
            obj.fracmesh = fracmesh;

            % Number of rock elements
            obj.num_rock_el = size(rockmesh.element.vertices,1);

            % Number of fracture elements and fracture indices
            [obj.num_frac_el, obj.frac_ind] = ...
                obj.define_frac_num(fracmesh);

            % Create fracture cell numbering dictionary
            [obj.frac_num, obj.frac_global, obj.frac_local] = ...
                obj.create_fracmesh_dictionary();

            % Initialize cross connections
            CrossConnections = obj.generate_empty_cross_connections();

            % Fill fracture-rock and fracture-fracture connections
            CrossConnections = ...
                obj.fillCrossConnections(CrossConnections, new_struct);

            % Initialize remaining cross connection properties
            obj.CrossConnections = fill_zeros(CrossConnections);

        end

        function [FracNum, FracGlobal, FracLocal] = ...
                create_fracmesh_dictionary(obj)

            % Global fracture cell numbering
            FracNum = repelem((1:length(obj.frac_ind))', obj.frac_ind);
            FracGlobal = (1:length(FracNum))';

            % Local fracture cell numbering
            FracLocal = zeros(sum(obj.num_frac_el),1);

            Accum = 0;
            Init = 1;

            % Loop through each fracture
            for Index = obj.num_frac_el'
                FracLocal((Init+Accum):(Index+Accum)) = (1:Index)';
                Accum = Accum + Index;
            end

        end

        function CrossConnections = ...
                fillCrossConnections(obj, CrossConnections, new_struct)

            % Offset for global fracture cell numbering
            LocSum = [0; obj.frac_ind];

            % Loop through EDFM fracture-rock connections
            for Index = 1:length(new_struct.fracrock_edfm)

                LocalEDFM = new_struct.fracrock_edfm{Index};

                % Sort connections by fracture cell
                [FracIndSorted, Ref] = sort(LocalEDFM.frac_ind);
                CISorted = LocalEDFM.CI(Ref);
                VolIndSorted = LocalEDFM.vol_ind(Ref);

                % Convert local fracture indices to global indices
                FracIndSorted = FracIndSorted + LocSum(Index);

                % Find unique fracture cells
                UniqueFracInd = unique(FracIndSorted)';

                % Loop through each fracture cell
                for LocalIndex = UniqueFracInd

                    FracElRef = FracIndSorted == LocalIndex;

                    % Add rock cell connections
                    CrossConnections(LocalIndex).Cells = ...
                        [CrossConnections(LocalIndex).Cells; ...
                         VolIndSorted(FracElRef)];

                    % Add connection indices
                    CrossConnections(LocalIndex).CI = ...
                        [CrossConnections(LocalIndex).CI; ...
                         CISorted(FracElRef)];

                end

            end

            if new_struct.pedfm_on

                % Loop through pEDFM fracture-rock connections
                for Index = 1:length(new_struct.fracrock_pedfm)

                    LocalPEDFM = new_struct.fracrock_pedfm{Index};

                    % Sort connections by fracture cell
                    [FracIndSorted, Ref] = ...
                        sort(LocalPEDFM.frac_ind);
                    CISorted = LocalPEDFM.CI(Ref);
                    VolIndSorted = LocalPEDFM.vol_ind(Ref);

                    % Convert local fracture indices to global indices
                    FracIndSorted = FracIndSorted + LocSum(Index);

                    % Loop through each fracture cell
                    for LocalIndex = unique(FracIndSorted)'

                        FracEl = LocalIndex;
                        Ref = FracIndSorted == FracEl;

                        % Add rock cell connections
                        CrossConnections(FracEl).Cells = ...
                            [CrossConnections(FracEl).Cells; ...
                             VolIndSorted(Ref)];

                        % Add connection indices
                        CrossConnections(FracEl).CI = ...
                            [CrossConnections(FracEl).CI; ...
                             CISorted(Ref)];

                        % Combine duplicate rock cell connections
                        [EDFMPEDFMCells,~,GroupsCell] = ...
                            unique(CrossConnections(FracEl).Cells);

                        CrossConnections(FracEl).CI = ...
                            accumarray(GroupsCell, ...
                                       CrossConnections(FracEl).CI);

                        CrossConnections(FracEl).Cells = ...
                            EDFMPEDFMCells;

                    end

                end

            end

            % Number of fractures
            nFrac = size(new_struct.fracfrac,1);

            % Loop through fracture-fracture connections
            for ii = 1:nFrac

                for jj = 1:nFrac

                    if ~isempty(new_struct.fracfrac{ii,jj})

                        Physical = new_struct.fracfrac{ii,jj};

                        % Convert local fracture indices to global indices
                        Ind1 = Physical.ind1;
                        CrossInd1 = Ind1 + LocSum(ii);

                        Ind2 = Physical.ind2;
                        CrossInd2 = Ind2 + LocSum(jj);

                        % Sort fracture-fracture connections
                        [~,Ref] = sortrows([CrossInd1 CrossInd2]);
                        CrossInd1 = CrossInd1(Ref);
                        CrossInd2 = CrossInd2(Ref);
                        CI = Physical.CI(Ref);

                        % Add connections from fracture 1 to fracture 2
                        for LocalIndex = 1:length(CrossInd1)

                            FracEl = CrossInd1(LocalIndex);
                            FracEl2 = CrossInd2(LocalIndex) + ...
                                      obj.num_rock_el;

                            CrossConnections(FracEl).Cells = ...
                                [CrossConnections(FracEl).Cells; ...
                                 FracEl2];

                            CrossConnections(FracEl).CI = ...
                                [CrossConnections(FracEl).CI; ...
                                 CI(LocalIndex)];

                        end

                        % Add connections from fracture 2 to fracture 1
                        for LocalIndex = 1:length(CrossInd1)

                            FracEl = CrossInd2(LocalIndex);
                            FracEl2 = CrossInd1(LocalIndex) + ...
                                      obj.num_rock_el;

                            CrossConnections(FracEl).Cells = ...
                                [CrossConnections(FracEl).Cells; ...
                                 FracEl2];

                            CrossConnections(FracEl).CI = ...
                                [CrossConnections(FracEl).CI; ...
                                 CI(LocalIndex)];

                        end

                        % Combine duplicate connections
                        TouchedFrac = ...
                            unique([CrossInd1(:); CrossInd2(:)]);

                        for t = 1:numel(TouchedFrac)

                            Fr = TouchedFrac(t);

                            if isempty(CrossConnections(Fr).Cells)
                                continue
                            end

                            Cells = CrossConnections(Fr).Cells(:);
                            CI = CrossConnections(Fr).CI(:);

                            [CellsUnique,~,Groups] = unique(Cells);
                            CIUnique = accumarray(Groups,CI);

                            [CellsUnique,SortIndex] = sort(CellsUnique);
                            CIUnique = CIUnique(SortIndex);

                            CrossConnections(Fr).Cells = CellsUnique;
                            CrossConnections(Fr).CI = CIUnique;

                        end

                    end

                end

            end

        end

        function CrossConnections = ...
                generate_empty_cross_connections(obj)

            % Initialize cross connection structure
            LCell = cell(1);

            CrossConnections = ...
                struct("Cells",LCell, ...
                       "CI",[], ...
                       "T_Geo",[], ...
                       "T_Geo_Cond",[], ...
                       "UpWind",[], ...
                       "U_Geo",[]);

            % Initialize one entry for each fracture cell
            CrossConnections(sum(obj.num_frac_el),1).CI = [];

        end

        function [NFracEl, FracElInd] = ...
                define_frac_num(~, fracmesh)

            % Number of fractures
            nFrac = length(fracmesh);

            % Initialize number of fracture cells
            NFracEl = zeros(nFrac,1);

            % Loop through each fracture
            for Index = 1:nFrac
                NFracEl(Index) = ...
                    length(fracmesh(Index).elements.area);
            end

            % Cumulative fracture cell indices
            FracElInd = cumsum(NFracEl);

        end

    end

end

function CrossConnections = fill_zeros(CrossConnections)

    % Loop through each fracture cell
    for Index = 1:length(CrossConnections)
    
        % Number of cross connections
        nCross = size(CrossConnections(Index).Cells,1);

        % Initialize geometric transmissibility
        CrossConnections(Index).T_Geo = ...
            zeros(nCross,1);
    
        % Initialize geometric conductive transmissibility
        CrossConnections(Index).T_Geo_Cond = ...
            zeros(nCross,1);
    
        % Initialize upwind values
        CrossConnections(Index).UpWind = ...
            zeros(nCross,2);
    
        % Initialize geometric velocity
        CrossConnections(Index).U_Geo = ...
            zeros(nCross,2);
    
    end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %DARSim Reservoir Simulator
% %Author: Artur Castiel;
% %TU Delft
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classdef translate_newDarsim_to_oldDarsim
%     properties
%         num_frac_el
%         num_rock_el
%         frac_ind
%         fracmesh
%         frac_num
%         frac_global
%         frac_local
%         CrossConnections
%     end
%     methods
%         function obj = translate_newDarsim_to_oldDarsim(new_struct, ...
%                 rockmesh, fracmesh)
%             obj.fracmesh = fracmesh;
%             obj.num_rock_el = size(rockmesh.element.vertices,1);
%             [obj.num_frac_el, obj.frac_ind] = ...
%                                         obj.define_frac_num(fracmesh);
%             [obj.frac_num,obj.frac_global,obj.frac_local] = ....
%                                 obj.create_fracmesh_dictionary();
% 
%             CrossConnections = generate_empty_cross_connections(obj);
%             [CrossConnections] = obj.fillCrossConnections(...
%                                        CrossConnections, new_struct);
%             obj.CrossConnections = fill_zeros(CrossConnections);
% 
%         end
% 
% 
%         function [frac_num,frac_global,frac_local] = ...
%                                 create_fracmesh_dictionary(obj)
%             frac_num = repelem([1:length(obj.frac_ind)]', obj.frac_ind);
%             frac_global = [1:length(frac_num)]';
%             frac_local = zeros(sum(obj.num_frac_el),1);
%             acum = 0;
%             init= 1 ;
%             for index = obj.num_frac_el'
%                 frac_local((init+acum):(index+acum)) = [1:index]';
%                 acum = acum + index;
%             end            
%         end
%         function CrossConnections = ....
%                 fillCrossConnections(obj,CrossConnections, new_struct)
%             %load Discretizations.mat
%             loc_sum = [0; obj.frac_ind];
%             for index = 1:length(new_struct.fracrock_edfm)
%                 local_edfm = new_struct.fracrock_edfm{index};
%                 [frac_ind_sorted,ref] = sort(local_edfm.frac_ind);
%                 CI_sorted = local_edfm.CI(ref);
%                 vol_ind_sorted = local_edfm.vol_ind(ref);
%                 %frac_ind_sorted2 = frac_ind_sorted;
%                 frac_ind_sorted = frac_ind_sorted + loc_sum(index);
%                 uniq_frac_ind = unique(frac_ind_sorted)';
% 
%                 for local_index = uniq_frac_ind %unique(frac_ind_sorted)'
%                     % if index == 1
%                     %     1
%                     % else
%                     %     1
%                     % end
%                     frac_el_ref = frac_ind_sorted == local_index;
%                     %frac_el = frac_ind_sorted(frac_el_ref);
% 
%                     CrossConnections(local_index).Cells = [...
%                         CrossConnections(local_index).Cells; ...
%                         vol_ind_sorted(frac_el_ref)];
%                     CrossConnections(local_index).CI = [...
%                         CrossConnections(local_index).CI;...
%                         CI_sorted(frac_el_ref)];
% 
%                     % CrossConnections(frac_el).Cells = [...
%                     %     CrossConnections(frac_el).Cells; ...
%                     %     vol_ind_sorted(local_index)];
%                     % CrossConnections(frac_el).CI = [...
%                     %     CrossConnections(frac_el).CI;...
%                     %     CI_sorted(local_index)];
%                 end
%             end
% 
%             if new_struct.pedfm_on
%                 for index = 1: length(new_struct.fracrock_pedfm)
%                     local_pedfm = new_struct.fracrock_pedfm{index};
%                     [frac_ind_sorted,ref] = sort(local_pedfm.frac_ind);
%                     CI_sorted = local_pedfm.CI(ref);
%                     vol_ind_sorted = local_pedfm.vol_ind(ref);
%                     frac_ind_sorted = frac_ind_sorted + loc_sum(index);
%                     for local_index = unique(frac_ind_sorted)'
%                         frac_el = local_index; %frac_ind_sorted(local_index);
%                         ref = frac_ind_sorted == frac_el;
%                         CrossConnections(frac_el).Cells = [...
%                             CrossConnections(frac_el).Cells; ...
%                             vol_ind_sorted(ref)];
%                         CrossConnections(frac_el).CI = [...
%                             CrossConnections(frac_el).CI;...
%                             CI_sorted(ref)];
%                         [edfm_pedfm_Cells,~,groupsCell] = unique(...
%                             CrossConnections(frac_el).Cells);
%                         CrossConnections(frac_el).CI = ...
%                             accumarray(groupsCell,CrossConnections(frac_el).CI);
%                         CrossConnections(frac_el).Cells = edfm_pedfm_Cells;
%                     end
%                 end
%             end
%             nfrac = size(new_struct.fracfrac,1);
%             for ii = 1:nfrac
%                 for jj = 1:nfrac
%                     if ~isempty(new_struct.fracfrac{ii,jj})
%                         physical = new_struct.fracfrac{ii,jj};
%                         ind1 = physical.ind1;
%                         cross_ind1 = ind1 + loc_sum(ii);
%                         ind2 = physical.ind2;
%                         cross_ind2 = ind2 + loc_sum(jj);
%                         [~,ref] = sortrows([cross_ind1 cross_ind2]);
%                         cross_ind1 = cross_ind1(ref);
%                         cross_ind2 = cross_ind2(ref);
%                         CI = physical.CI(ref);
%                         for local_index = 1:length(cross_ind1)
%                             frac_el = cross_ind1(local_index);
%                             frac_el2 = cross_ind2(local_index) + ...
%                                 obj.num_rock_el;
%                             CrossConnections(frac_el).Cells = [...
%                                 CrossConnections(frac_el).Cells; frac_el2];
%                             CrossConnections(frac_el).CI = [...
%                                 CrossConnections(frac_el).CI;...
%                                 CI(local_index)];
%                         end
%                         for local_index = 1:length(cross_ind1)
%                             frac_el = cross_ind2(local_index);
%                             frac_el2 = cross_ind1(local_index) + ...
%                                 obj.num_rock_el;
%                             CrossConnections(frac_el).Cells = [...
%                                 CrossConnections(frac_el).Cells; frac_el2];
%                             CrossConnections(frac_el).CI = [...
%                                 CrossConnections(frac_el).CI;...
%                                 CI(local_index)];
%                         end
% 
%                         touched_frac = unique([cross_ind1(:); cross_ind2(:)]);
%                         for t = 1:numel(touched_frac)
%                             fr = touched_frac(t);
%                             if isempty(CrossConnections(fr).Cells), continue; end
%                             cells = CrossConnections(fr).Cells(:);
%                             ci    = CrossConnections(fr).CI(:);
%                             [cells_u, ~, grp] = unique(cells);
%                             ci_u = accumarray(grp, ci);
%                             [cells_u, p] = sort(cells_u);
%                             ci_u = ci_u(p);
% 
%                             CrossConnections(fr).Cells = cells_u;
%                             CrossConnections(fr).CI    = ci_u;
%                         end
%                     end
%                 end
%             end
%         end
%         function [CrossConnections] = generate_empty_cross_connections(obj)
%             lcell = cell(1);
%             CrossConnections = struct("Cells",lcell, "CI", [], ...
%                 "T_Geo", [], "T_Geo_Cond",[], "UpWind",[], "U_Geo", []);
%             CrossConnections(sum(obj.num_frac_el),1).CI = [];
%         end
% 
%         function [frac_num_el,frac_el_ind] = define_frac_num(~, fracmesh)
%             num_frac =length(fracmesh);
%             frac_num_el = zeros(num_frac,1);
% 
%             for index = 1:num_frac
%                 frac_num_el(index) = length(fracmesh(index).elements.area);                
%             end
%             frac_el_ind = cumsum(frac_num_el);
%         end
% 
%     end
% end
% 
% 
% function [CrossConnections] = fill_zeros(CrossConnections)
% for index = 1:length(CrossConnections)
%     num_cross = size(CrossConnections(index).Cells,1);
%     CrossConnections(index).T_Geo = zeros(num_cross,1);
%     CrossConnections(index).T_Geo_Cond = ...
%         zeros(num_cross,1);
%     CrossConnections(index).UpWind = ...
%         zeros(num_cross,2);
%     CrossConnections(index).U_Geo = ...
%         zeros(num_cross,2);
% end
% end
