%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Finding EDFM and pEDFM connections and intersection areas. 
% Intersection areas are calculated geometrically.
%
% Input:
%   FracMesh:       struct
%                   Contains relevant frac mesh information
%   RockMesh:       struct
%                   Contains relevant rock mesh information
%   FracGeometry:   struct
%                   Contains relevant information about fracture properties
%   Settings:       struct
%                   Contains relevant information about simulation settings
%
% Output:
%   EDFMFracRock:   struct
%                   Contains three cell fields:
%                   RockID - rock cell index of EDFM connection
%                   FracID - fracture cell index of EDFM connection
%                   Area - area of EDFM connection
%   PEDFMFracRock:  struct
%                   Contains four cell field: 
%                   RockID - rock cell index of pEDFM connection
%                   FracID - fracture cell index of pEDFM connection
%                   Area - area of pEDFM connection
%                   AlphaID - face index with altered alpha values in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EDFMFracRock, ...
          PEDFMFracRock ] = ...
          FracRockIntersection(FracMesh, ...
                               RockMesh, ...
                               FracGeometry, ...
                               Settings)

% Number of fractures in simulation
nFrac = max(size(FracMesh));

% Initialize output variables
EDFMFracRock = struct();
RockCellEDFMOutput = cell(nFrac,1);
FracCellEDFMOutput = cell(nFrac,1);
AreaEDFMOutput = cell(nFrac,1);
PEDFMFracRock = struct();
RockCellPEDFMOutput = cell(nFrac,1);
FracCellPEDFMOutput = cell(nFrac,1);
AreaPEDFMOutput = cell(nFrac,1);
AlphaFaceOutput = cell(nFrac,1);

% Loop over each fracture individually
for iFrac = 1:nFrac

    % Find intesecting rock cells and lengths of intersections (EDFM)
    RockIntersections = find_edfm_frac_rock(RockMesh, ... 
                                            FracMesh,...
                                            iFrac, ...
                                            FracGeometry(iFrac));

    % Place in output Cell objects
    RockCellEDFMOutput{iFrac,1} = RockIntersections.RockCellID;
    FracCellEDFMOutput{iFrac,1} = RockIntersections.FracCellID;
    AreaEDFMOutput{iFrac,1} = RockIntersections.AreaEDFM;
    
    % Find pEDFM conenctions
    if Settings.pedfm
        EDFMconnections = ...
            [RockIntersections.RockCellID, RockIntersections.FracCellID];
        [PEDFMconnections, AreaPEDFM, AlphaFaceID]= ...
            find_pdefm_frac_rock(RockIntersections, ...
                                 RockMesh, ...
                                 FracGeometry(iFrac));
        RockCellPEDFMOutput{iFrac,1} = PEDFMconnections(:,1); 
        FracCellPEDFMOutput{iFrac,1} = PEDFMconnections(:,2);
        AreaPEDFMOutput{iFrac,1} = AreaPEDFM;  
        AlphaFaceOutput{iFrac,1} = AlphaFaceID;
    end 

end

% Save Output
EDFMFracRock.RockID = RockCellEDFMOutput;
EDFMFracRock.FracID = FracCellEDFMOutput;
EDFMFracRock.Area = AreaEDFMOutput;

if Settings.pedfm
    % Group alpha values
    AlphaFaceOutput = cell2mat(AlphaFaceOutput);
    [AlphaFaceOutput,~,~] = unique(AlphaFaceOutput);

    PEDFMFracRock.RockID = RockCellPEDFMOutput;
    PEDFMFracRock.FracID = FracCellPEDFMOutput;
    PEDFMFracRock.Area = AreaPEDFMOutput;
    PEDFMFracRock.AlphaID = AlphaFaceOutput;
end

end