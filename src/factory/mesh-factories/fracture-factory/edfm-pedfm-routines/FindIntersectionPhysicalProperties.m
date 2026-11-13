%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% Calculating the physical properties of EDFM and pEDFM connections between
% fracture-rock cells and fracture-fracture cells.
%
% Input:
%   RockMesh:           struct
%                       Contains relevant rock mesh information
%   FracMesh:           struct
%                       Contains relevant fracture mesh information
%   FracGeometry:       struct
%                       Contains relevant fracture geometry and properties
%   LFracFracEDFM:      cell
%                       Contains length of frac-frac intersections EDFM
%   LFracFracPEDFM:     cell
%                       Contains length of frac-frac intersections pEDFM
%   AFracRockEDFM:      struct
%                       Contains area of frac-rock intersections EDFM
%   AFrackRockPEDFM:    struct
%                       Contains area of frac-rock intersections pEDFM
%   DFracFrac:          cell
%                       Contains distances between frac-frac connections
%                       cells
%   DFracRock:          cell
%                       Contains distances between frac-rock connections
%   Settings:           struct
%                       Contains relevant simulation settings
%
% Output:
%   OutputPhysProp:     struct
%                       Contains physical properties for intra-fracture,
%                       fracture-fracture, EDFM, and pEDFM connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OutputPhysProp] = ...
         FindIntersectionPhysicalProperties(RockMesh, ...
                                            FracMesh, ...
                                            FracGeometry, ...
                                            LFracFracEDFM, ...
                                            LFracFracPEDFM, ...
                                            AFracRockEDFM, ...
                                            AFrackRockPEDFM, ...
                                            DFracFrac, ...
                                            DFracRock, ...
                                            Settings)

% Initialize output structure
OutputPhysProp = struct();

% Number of fractures
nFrac = size(FracMesh,1);

% Calculate fracture-fracture EDFM connection properties
FracFracEDFMProp = cell(nFrac,nFrac);

for iFrac1 = 1:nFrac

    for iFrac2 = iFrac1+1:nFrac

        if ~isempty(LFracFracEDFM.Length{iFrac1,iFrac2})

            % Initialize physical properties
            PhysProp = struct();

            % Fracture cells connected through fracture intersection
            LocalFrac1ID = LFracFracEDFM.Frac1ID{iFrac1,iFrac2};
            LocalFrac2ID = LFracFracEDFM.Frac2ID{iFrac1,iFrac2};

            % Length of fracture-fracture intersection within each
            % connected pair of fracture cells
            LocalLengthEDFM = LFracFracEDFM.Length{iFrac1,iFrac2};

            % Average distance from fracture cells to intersection line
            DistInfoFrac1CellID = DFracFrac{iFrac1}{iFrac2,1};
            DistInfoFrac2CellID = DFracFrac{iFrac2}{iFrac1,1};
            DistInfoFrac1Dist = DFracFrac{iFrac1}{iFrac2,2};
            DistInfoFrac2Dist = DFracFrac{iFrac2}{iFrac1,2};

            % Map average distances to each fracture-fracture connection
            [~,idx] = ismember(LocalFrac1ID,DistInfoFrac1CellID);
            LocalFrac1Dist = DistInfoFrac1Dist(idx);
            [~,idx] = ismember(LocalFrac2ID,DistInfoFrac2CellID);
            LocalFrac2Dist = DistInfoFrac2Dist(idx);

            % Fracture permeabilities
            K1 = FracGeometry(iFrac1).Permeability;
            K2 = FracGeometry(iFrac2).Permeability;

            % Harmonic average permeability between fractures
            PhysProp.K = (2*K1*K2)./(K1+K2);

            % Store EDFM connections
            PhysProp.Frac1ID = LocalFrac1ID;
            PhysProp.Frac2ID = LocalFrac2ID;

            % Connectivity index from each fracture to intersection line
            CI1 = 2*LocalLengthEDFM./LocalFrac1Dist;
            CI2 = 2*LocalLengthEDFM./LocalFrac2Dist;

            % Harmonic average fracture-fracture connectivity index
            CIHarm = 2*(CI1.*CI2)./(CI1+CI2);
            CIHarm(isnan(CIHarm)) = 0;
            PhysProp.CI = CIHarm;

            FracFracEDFMProp{iFrac1,iFrac2} = PhysProp;

        end

    end

end

OutputPhysProp.FracFracEDFM = FracFracEDFMProp;

% Calculate fracture-rock EDFM connection properties
FracRockEDFMProp = cell(nFrac,1);

for iFrac = 1:nFrac

    % Initialize physical properties
    PhysProp = struct();

    % Rock and fracture cell indices of EDFM connections
    LocalRockID = AFracRockEDFM.RockID{iFrac};
    LocalFracID = AFracRockEDFM.FracID{iFrac};

    % Fracture-rock intersection area
    LocalAreaEDFM = AFracRockEDFM.Area{iFrac};

    % Average distance between fracture and rock cells
    DistInfoRockCellID = DFracRock{iFrac,1};
    DistInfoRockDist = DFracRock{iFrac,2};

    % Map average distances to each fracture-rock connection
    [~,idx] = ismember(LocalRockID, DistInfoRockCellID);
    LocalRockDist = DistInfoRockDist(idx);

    % EDFM connectivity index
    CI = LocalAreaEDFM./LocalRockDist;

    % Rock and fracture permeability
    K1 = RockMesh.props.perm(LocalRockID);
    K2 = FracGeometry(iFrac).Permeability;
    K2 = repmat(K2,size(K1,1),1);

    % Characteristic rock cell length and fracture aperture
    L = RockMesh.element.volume(LocalRockID).^(1/3);
    W = FracGeometry(iFrac).Aperture;

    % Store EDFM connection properties
    PhysProp.RockID = LocalRockID;
    PhysProp.FracID = LocalFracID;

    % Effective permeability between rock and fracture cells
    PhysProp.K = (L+W).*(K1.*K2)./(W.*K1+L.*K2);
    PhysProp.CI = CI;

    FracRockEDFMProp{iFrac} = PhysProp;

end

OutputPhysProp.FracRockEDFM = FracRockEDFMProp;

% Initialize pEDFM properties
OutputPhysProp.PEDFM_on = logical(Settings.pedfm);

% Calculate pEDFM connection properties
if OutputPhysProp.PEDFM_on

    % Calculate the frac-frac pEDFM connection properties
    FracFracPEDFMProp = cell(nFrac,1);

    for iFrac1 = 1:nFrac

        for iFrac2 = iFrac1+1:nFrac

            if ~isempty(LFracFracPEDFM.Length{iFrac1,iFrac2})

                % Initialize physical properties
                PhysProp = struct();
                
                % Fracture cell indices of pEDFM connections
                LocalFrac1ID = LFracFracPEDFM.Frac1ID{iFrac1,iFrac2};
                LocalFrac2ID = LFracFracPEDFM.Frac2ID{iFrac1,iFrac2};
    
                % Length of frac-frac intersection for pEDFM connection
                LocalLengthPEDFM = LFracFracPEDFM.Length{iFrac1,iFrac2};

                % Average distance from fracture cells to intersection line
                DistInfoFrac1CellID = DFracFrac{iFrac1}{iFrac2,1};
                DistInfoFrac2CellID = DFracFrac{iFrac2}{iFrac1,1};
                DistInfoFrac1Dist = DFracFrac{iFrac1}{iFrac2,2};
                DistInfoFrac2Dist = DFracFrac{iFrac2}{iFrac1,2};
    
                % Map average distances to each pEDFM connection
                [~,idx] = ismember(LocalFrac1ID,DistInfoFrac1CellID);
                LocalFrac1Dist = DistInfoFrac1Dist(idx);
                [~,idx] = ismember(LocalFrac2ID,DistInfoFrac2CellID);
                LocalFrac2Dist = DistInfoFrac2Dist(idx);

                % Fracture permeabilities
                K1 = FracGeometry(iFrac1).Permeability;
                K2 = FracGeometry(iFrac2).Permeability;
    
                % Harmonic average permeability between fractures
                PhysProp.K = (2*K1*K2)./(K1+K2);

                % Store pEDFM connections
                PhysProp.Frac1ID = LocalFrac1ID;
                PhysProp.Frac2ID = LocalFrac2ID;
    
                % Connectivity index from each frac to intersection line
                CI1 = LocalLengthPEDFM./LocalFrac1Dist;
                CI2 = LocalLengthPEDFM./LocalFrac2Dist;
    
                % Harmonic average fracture-fracture connectivity index
                CIHarm = 2*(CI1.*CI2)./(CI1+CI2);
                CIHarm(isnan(CIHarm)) = 0;
                PhysProp.CI = CIHarm;

                FracFracPEDFMProp{iFrac1,iFrac2} = PhysProp;

            end

        end

    end

    OutputPhysProp.FracFracPEDFM = FracFracPEDFMProp;
    
    % Calculate the frac-rock pEDFM connection properties
    FracRockPEDFMProp = cell(nFrac,1);

    for iFrac = 1:nFrac

        % Initialize physical properties
        PhysProp = struct();

        % Rock and fracture cell indices of pEDFM connections
        LocalRockID = AFrackRockPEDFM.RockID{iFrac};
        LocalFracID = AFrackRockPEDFM.FracID{iFrac};

        % Projected fracture-rock intersection area
        LocalAreaEDFM = AFrackRockPEDFM.Area{iFrac};

        % Average distance between fracture and rock cells
        DistInfoRockCellID = DFracRock{iFrac,1};
        DistInfoRockDist = DFracRock{iFrac,2};

        % Map average distances to each fracture-rock connection
        [~,idx] = ismember(LocalRockID, DistInfoRockCellID);
        LocalRockDist = DistInfoRockDist(idx);

        % pEDFM connectivity index
        CI = LocalAreaEDFM./LocalRockDist;

        % Rock and fracture permeability
        K1 = RockMesh.props.perm(LocalRockID);
        K2 = FracGeometry(iFrac).Permeability;
        K2 = repmat(K2,size(K1,1),1);

        % Characteristic rock cell length and fracture aperture
        L = RockMesh.element.volume(LocalRockID).^(1/3);
        W = FracGeometry(iFrac).Aperture;
        
        % Store pEDFM connection properties
        PhysProp.RockID = LocalRockID;
        PhysProp.FracID = LocalFracID;

        % Effective permeability between rock and fracture cells
        PhysProp.K = (L+W).*(K1.*K2)./(W.*K1+L.*K2);
        PhysProp.CI = CI;

        FracRockPEDFMProp{iFrac} = PhysProp;

    end

    OutputPhysProp.FracRockPEDFM = FracRockPEDFMProp;

end

% Calculate pEDFM alpha values
OutputPhysProp.RockAlpha = zeros(RockMesh.face.inner.num,1);
OutputPhysProp.FracAlpha = cell(nFrac,1);

if OutputPhysProp.PEDFM_on

    % Frac Mesh Alpha
    if nFrac > 1
        for iFrac = 1:nFrac
    
            % Number of faces in frac mesh
            nFaces = FracMesh(iFrac).faces.inner.num + ...
                     FracMesh(iFrac).faces.bound.num;
    
            % Setting alpha to 0 (no projection) or 1 (projection)
            % In this method, we are assuming alpha = 1 for all projection faces in
            % fracture projection path
            AlphaGeom = zeros(nFaces,1);
            AlphaGeom(LFracFracPEDFM.AlphaID{iFrac}) = 1;
    
            OutputPhysProp.FracAlpha{iFrac} = AlphaGeom;
    
        end
    end
    
    % Rock Mesh Alpha
    % Number of faces in rock mesh
    nFaces = RockMesh.face.inner.num + ...
             RockMesh.face.bound.num;

    % Setting alpha to 0 (no projection) or 1 (projection)
    % In this method, we are assuming alpha = 1 for all projection faces in
    % fracture projection path
    AlphaGeom = zeros(nFaces,1);
    AlphaGeom(AFrackRockPEDFM.AlphaID(:)) = 1;

    % Store alpha values
    OutputPhysProp.RockAlpha = AlphaGeom;

end

end