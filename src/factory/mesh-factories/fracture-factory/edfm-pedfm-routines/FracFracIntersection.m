%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Finding EDFM and pEDFM connections and intersection lengths. 
% Intersection lengths are calculated geometrically.
%
% Input:
%   FracMesh:           struct
%                       Contains relevant frac mesh information
%   GeoIntersection:    struct
%                       Contains relevant info on fracture intersections
%   Settings:       struct
%                   Contains relevant information about simulation settings
%
% Output:
%   EDFMFracFrac:   struct; with three fields:
%                   "Frac1ID" (f1 cell ID)
%                   "Frac2ID" (f2 cell ID)
%                   "Length" (Length of intersection line)
%                   Each field contains cell (NumFrac-by-NumFrac)
%                   Row represents f1, column represents f2
%                   Each entry contains a list of values for specific field
%                   e.g. EDFMFracFrac.Frac1ID{3,2} gives a list of frac
%                   cell IDs from frac 3 that intersect with frac 2.
%   PEDFMFracFrac:  struct; with four fields:
%                   "Frac1ID" (f1 cell ID)
%                   "Frac2ID" (f2 cell ID)
%                   "Length" (Length of intersection line)
%                   "AlphaID" (Edge/Face IDs in frace mesh for alpha proj.)
%                   First three fields contains cell (NumFrac-by-NumFrac)
%                   Row represents f1, column represents f2
%                   Each entry contains a list of values for specific field
%                   e.g. PEDFMFracFrac.Frac1ID{3,2} gives a list of frac
%                   cell IDs from frac 3 with PEDFM connections to frac 2. 
%                   Last field (AlphaID) contains cell (NumFrac-by-1)
%                   Row represents frac ID
%                   Each entry contains a list of face IDs in local frac
%                   mesh for alpha projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EDFMFracFrac, ...
          PEDFMFracFrac] = ...
          FracFracIntersection(FracMesh, ...
                               GeoIntersection, ...
                               Settings)

% Number of fractures in simulation
nFrac = max(size(FracMesh));

% Initialize output variables
EDFMFracFrac = struct();
Frac1CellEDFMOutput = cell(nFrac,nFrac);
Frac2CellEDFMOutput = cell(nFrac,nFrac);
LengthEDFMOutput = cell(nFrac,nFrac);
PEDFMFracFrac = struct();
Frac1CellPEDFMOutput = cell(nFrac,nFrac);
Frac2CellPEDFMOutput = cell(nFrac,nFrac);
LengthPEDFMOutput = cell(nFrac,nFrac);
AlphaFaceOutput = cell(nFrac,1);

% Loop over each fracture individually
for iFrac1 = 1:nFrac

    % Find which fractures this fracture intersects with
    IntersectsWithThisFrac = find(GeoIntersection.intersectWith(:,iFrac1));
    
    % Loop through other fractures with index greater than this one
    for iFrac2 = 1:nFrac
       
        % Check if this fracture intersects with iFrac1
        if ismember(iFrac2,IntersectsWithThisFrac)

            % Find the intersection line
            IntLineEndPoints = cell2mat(...
                GeoIntersection.get_boundary_line_intersection(iFrac1, ...
                                                               iFrac2));

            % Find intesecting fracture cells and lengths of intersections
            FractureIntersections = ...
                find_edfm_frac_frac(FracMesh, ...
                                    iFrac1, ...
                                    iFrac2, ...
                                    IntLineEndPoints);

            % Fill in output
            Frac1CellEDFMOutput{iFrac1,iFrac2} = ...
                FractureIntersections.CellID1;
            Frac2CellEDFMOutput{iFrac1,iFrac2} = ...
                FractureIntersections.CellID2;
            LengthEDFMOutput{iFrac1,iFrac2} = ...
                FractureIntersections.LengthEDFM;
            
            % Find pEDFM conenctions
            if Settings.pedfm
                [PEDFMCellIDs, PEDFMLength, AlphaFaceIDs] = ...
                    find_pdefm_frac_frac(FractureIntersections, ...
                                         FracMesh, ...
                                         iFrac1, ...
                                         IntLineEndPoints);
                if ~isempty(PEDFMCellIDs)
                    Frac1CellPEDFMOutput{iFrac1,iFrac2} = ...
                        PEDFMCellIDs(:,1);
                    Frac2CellPEDFMOutput{iFrac1,iFrac2} = ...
                        PEDFMCellIDs(:,2);
                    LengthPEDFMOutput{iFrac1,iFrac2} = ...
                        PEDFMLength;
                end
                if ~isempty(AlphaFaceIDs)
                    AlphaFaceOutput{iFrac1} = AlphaFaceIDs;
                end
            end
            
        end

    end

end

% Clean up output
% Fill in 'NA' for f1 => f2 for EDFM output
for iFrac1 = 1:nFrac
    for iFrac2 = 1:nFrac
        if iFrac1 >= iFrac2
            Frac1CellEDFMOutput{iFrac1,iFrac2} = 'NA';
            Frac2CellEDFMOutput{iFrac1,iFrac2} = 'NA';
            LengthEDFMOutput{iFrac1,iFrac2} = 'NA';
        end
    end
end

% Collect output for certain f1, f2 pairs for pEDFM output
% Fill in 'NA' for f1 => f2 for pEDFM output
for iFrac1 = 1:nFrac
    for iFrac2 = 1:nFrac
        if iFrac1 == iFrac2
            Frac1CellPEDFMOutput{iFrac1,iFrac2} = 'NA';
            Frac2CellPEDFMOutput{iFrac1,iFrac2} = 'NA';
            LengthPEDFMOutput{iFrac1,iFrac2} = 'NA';
        elseif iFrac1 > iFrac2
            Frac1CellPEDFMOutput{iFrac2,iFrac1} = ...
                [Frac1CellPEDFMOutput{iFrac2,iFrac1}; ...
                 Frac2CellPEDFMOutput{iFrac1,iFrac2}];
            Frac2CellPEDFMOutput{iFrac2,iFrac1} = ...
                [Frac2CellPEDFMOutput{iFrac2,iFrac1}; ...
                 Frac1CellPEDFMOutput{iFrac1,iFrac2}];
            LengthPEDFMOutput{iFrac2,iFrac1} = ...
                [LengthPEDFMOutput{iFrac2,iFrac1}; ...
                 LengthPEDFMOutput{iFrac1,iFrac2}];
            Frac1CellPEDFMOutput{iFrac1,iFrac2} = 'NA';
            Frac2CellPEDFMOutput{iFrac1,iFrac2} = 'NA';
            LengthPEDFMOutput{iFrac1,iFrac2} = 'NA';
        end
    end
end


% Save Output
EDFMFracFrac.Frac1ID = Frac1CellEDFMOutput;
EDFMFracFrac.Frac2ID = Frac2CellEDFMOutput;
EDFMFracFrac.Length = LengthEDFMOutput;

if Settings.pedfm
    PEDFMFracFrac.Frac1ID = Frac1CellPEDFMOutput;
    PEDFMFracFrac.Frac2ID = Frac2CellPEDFMOutput;
    PEDFMFracFrac.Length = LengthPEDFMOutput;
    PEDFMFracFrac.AlphaID = AlphaFaceOutput;
end

end

