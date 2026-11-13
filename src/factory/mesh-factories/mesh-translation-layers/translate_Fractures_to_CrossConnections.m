%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% Translating new DARSim fracture connectivity data to the old DARSim
% cross connection format.
%
% Input:
%   FracCore:       struct
%                   Contains fracture-rock and fracture-fracture
%                   connectivity information from the new DARSim format
%   RockMesh:       struct
%                   Contains relevant rock mesh information
%   FracMesh:       struct
%                   Contains relevant fracture mesh information
%
% Output:
%   CrossConnections: struct
%                     Contains cross connection information for each
%                     fracture cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef translate_Fractures_to_CrossConnections

    properties
        FracMesh
        nFrac
        nRockElem
        nFracElem
        CummFracElem
        FracIDByFracCell
        FracCellIDGlobal
        FracCellIDLocal
        CrossConnections
    end

    methods

        function obj = ...
                translate_Fractures_to_CrossConnections(FracCore, ...
                                                        RockMesh, ...
                                                        FracMesh)

            % Store fracture mesh
            obj.FracMesh = FracMesh;

            % Number of fractures
            obj.nFrac = max(size(FracMesh));

            % Number of rock elements
            obj.nRockElem = size(RockMesh.element.vertices,1);

            % Number of fracture elements and fracture indices
            [obj.nFracElem, obj.CummFracElem] = ...
                obj.define_frac_num(FracMesh);

            % Create fracture cell numbering dictionary
            [obj.FracIDByFracCell, ...
             obj.FracCellIDGlobal, ...
             obj.FracCellIDLocal] = obj.create_fracmesh_dictionary();

            % Initialize cross connections
            CrossConnections = obj.generate_empty_cross_connections();

            % Fill fracture-rock and fracture-fracture connections
            CrossConnections = ...
                obj.fillCrossConnections(CrossConnections, FracCore);

            % Initialize remaining cross connection properties
            obj.CrossConnections = fill_zeros(CrossConnections);

        end

        function [FracNum, FracGlobal, FracLocal] = ...
                create_fracmesh_dictionary(obj)

            % Global fracture cell numbering: In frac space (w/o rock mesh)
            FracNum = ...
                repelem((1:length(obj.CummFracElem))', obj.nFracElem);
            FracGlobal = (1:length(FracNum))';

            % Local fracture cell numbering
            FracLocal = zeros(sum(obj.nFracElem),1);

            Accum = 0;
            Init = 1;

            % Loop through each fracture
            for Index = obj.nFracElem'
                FracLocal((Init+Accum):(Index+Accum)) = (1:Index)';
                Accum = Accum + Index;
            end

        end

        function CrossConnections = ...
                fillCrossConnections(obj, CrossConnections, FracCore)

            % Offset for global fracture cell numbering
            % i.e. converting from local frac space to global frac space
            % This indexing still ignores the rock mesh
            LocSum = [0; obj.CummFracElem];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % EDFM Frac-Rock connections %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Loop through fractures
            for iFrac = 1:length(FracCore.FracRockEDFM)

                LocalEDFM = FracCore.FracRockEDFM{iFrac};

                % Sort connections by fracture cell
                [FracIndSorted, FracElRef] = sort(LocalEDFM.FracID);
                CISorted = LocalEDFM.CI(FracElRef);
                RockIndSorted = LocalEDFM.RockID(FracElRef);

                % Convert local frac elem ID to global frac space
                FracIndSorted = FracIndSorted + LocSum(iFrac);

                % Find unique fracture cells
                UniqueFracInd = unique(FracIndSorted)';

                % Loop through each fracture cell
                for LocalIndex = UniqueFracInd

                    FracElRef = FracIndSorted == LocalIndex;

                    % Add rock cell connections
                    CrossConnections(LocalIndex).Cells = ...
                        [CrossConnections(LocalIndex).Cells; ...
                         RockIndSorted(FracElRef)];

                    % Add connection indices
                    CrossConnections(LocalIndex).CI = ...
                        [CrossConnections(LocalIndex).CI; ...
                         CISorted(FracElRef)];

                end

            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PEDFM Frac-Rock connections %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if FracCore.PEDFM_on

                % Loop through fracs
                for iFrac = 1:length(FracCore.FracRockPEDFM)

                    LocalPEDFM = FracCore.FracRockPEDFM{iFrac};

                    % Sort connections by fracture cell
                    [FracIndSorted, FracElRef] = ...
                        sort(LocalPEDFM.FracID);
                    CISorted = LocalPEDFM.CI(FracElRef);
                    RockIndSorted = LocalPEDFM.RockID(FracElRef);

                    % Convert local frac elem ID to global frac space
                    FracIndSorted = FracIndSorted + LocSum(iFrac);
                    
                    % Find unique fracture cells
                    UniqueFracInd = unique(FracIndSorted)';

                    % Loop through each fracture cell
                    for LocalIndex = UniqueFracInd

                        FracElRef = FracIndSorted == LocalIndex;

                        % Add rock cell connections
                        CrossConnections(LocalIndex).Cells = ...
                            [CrossConnections(LocalIndex).Cells; ...
                             RockIndSorted(FracElRef)];

                        % Add connection indices
                        CrossConnections(LocalIndex).CI = ...
                            [CrossConnections(LocalIndex).CI; ...
                             CISorted(FracElRef)];

                    end

                end

            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % EDFM Frac-Frac connections %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Loop through fracture-fracture connections
            for iFrac1 = 1:obj.nFrac
                
                % Don't double count EDFM connections with fractures that
                % have already been read in with "for iFrac1 = 1:obj.nFrac"
                for iFrac2 = iFrac1+1:obj.nFrac

                    if ~isempty(FracCore.FracFracEDFM{iFrac1,iFrac2})

                        PhysProp = FracCore.FracFracEDFM{iFrac1,iFrac2};

                        % Convert local frac elem ID to global frac space
                        Frac1ID = PhysProp.Frac1ID;
                        Frac1ID = Frac1ID + LocSum(iFrac1);

                        Frac2ID = PhysProp.Frac2ID;
                        Frac2ID = Frac2ID + LocSum(iFrac2);

                        % Sort fracture-fracture connections
                        [~,FracElRef] = sortrows([Frac1ID Frac2ID]);
                        Frac1ID = Frac1ID(FracElRef);
                        Frac2ID = Frac2ID(FracElRef);
                        CI = PhysProp.CI(FracElRef);

                        % Add connections from fracture 1 to fracture 2
                        % We only need to add the EDFM connections from one
                        % frac to another once, so from fracture 1 to
                        % fracture 2 is sufficient. Don't add from fracture
                        % 2 to fracture 1, as this would double count EDFM
                        % connections. 
                        for LocalIndex = 1:length(Frac1ID)
                            
                            % Frac IDs in global fracture space
                            FracEl1 = Frac1ID(LocalIndex);
                            FracEl2 = Frac2ID(LocalIndex);
                            
                            % Frac IDs in global space, including rock mesh
                            FracEl2FullGlobal = FracEl2 + obj.nRockElem;

                            CrossConnections(FracEl1).Cells = ...
                                [CrossConnections(FracEl1).Cells; ...
                                 FracEl2FullGlobal];

                            CrossConnections(FracEl1).CI = ...
                                [CrossConnections(FracEl1).CI; ...
                                 CI(LocalIndex)];

                        end

                    end

                end

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PEDFM Frac-Frac connections %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Loop through fracture-fracture connections
            for iFrac1 = 1:obj.nFrac
                
                % Don't double count PEDFM connections with fractures that
                % have already been read in with "for iFrac1 = 1:obj.nFrac"
                for iFrac2 = iFrac1+1:obj.nFrac

                    if ~isempty(FracCore.FracFracPEDFM{iFrac1,iFrac2})

                        PhysProp = FracCore.FracFracPEDFM{iFrac1,iFrac2};

                        % Convert local frac elem ID to global frac space
                        Frac1ID = PhysProp.Frac1ID;
                        Frac1ID = Frac1ID + LocSum(iFrac1);

                        Frac2ID = PhysProp.Frac2ID;
                        Frac2ID = Frac2ID + LocSum(iFrac2);

                        % Sort fracture-fracture connections
                        [~,FracElRef] = sortrows([Frac1ID Frac2ID]);
                        Frac1ID = Frac1ID(FracElRef);
                        Frac2ID = Frac2ID(FracElRef);
                        CI = PhysProp.CI(FracElRef);

                        % Add connections from fracture 1 to fracture 2
                        % We only need to add the pEDFM connections from 
                        % one frac to another once, so from fracture 1 to
                        % fracture 2 is sufficient. Don't add from fracture
                        % 2 to fracture 1, as this would double count pEDFM
                        % connections. 
                        for LocalIndex = 1:length(Frac1ID)
                            
                            % Frac IDs in global fracture space
                            FracEl1 = Frac1ID(LocalIndex);
                            FracEl2 = Frac2ID(LocalIndex);
                            
                            % Frac IDs in global space, including rock mesh
                            FracEl2FullGlobal = FracEl2 + obj.nRockElem;

                            CrossConnections(FracEl1).Cells = ...
                                [CrossConnections(FracEl1).Cells; ...
                                 FracEl2FullGlobal];

                            CrossConnections(FracEl1).CI = ...
                                [CrossConnections(FracEl1).CI; ...
                                 CI(LocalIndex)];

                        end

                    end

                end

            end

            % Combine duplicate connections
            for LocalIndex = obj.FracCellIDGlobal'

                [EDFMPEDFMCells,~,GroupsCell] = ...
                    unique(CrossConnections(LocalIndex).Cells);
    
                CrossConnections(LocalIndex).CI = ...
                    accumarray(GroupsCell, ...
                               CrossConnections(LocalIndex).CI);
    
                CrossConnections(LocalIndex).Cells = ...
                    EDFMPEDFMCells;
            
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
            CrossConnections(sum(obj.nFracElem),1).CI = [];

        end

        function [nFracEl, CummFracElem] = ...
                define_frac_num(obj, fracmesh)

            % Initialize number of fracture cells
            nFracEl = zeros(obj.nFrac,1);

            % Loop through each fracture
            for Index = 1:obj.nFrac
                nFracEl(Index) = ...
                    length(fracmesh(Index).elements.area);
            end

            % Cumulative fracture cell indices
            CummFracElem = cumsum(nFracEl);

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
