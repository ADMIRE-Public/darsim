%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% Monte Carlo alternative to FracRockIntersection (EDFM_METHOD
% monte_carlo). Test points sampled on every fracture cell are assigned to
% the rock cells that contain them; hit ratios give the intersection areas
% (EDFM). pEDFM projection faces are selected as in find_pdefm_frac_rock
% (fracture plane crossing the centroid-centroid connector), applied to the
% sampled EDFM areas. Settings.pedfm_projection = "sampled" uses instead
% the nearest face per axis with the share of projected sample points.
%
% Input / Output: same contract as FracRockIntersection, so that
% find_distance_frac_rock(_MC) and FindIntersectionPhysicalProperties can
% consume the result unchanged.
%   EDFMFracRock.Area  is two-sided (2 x one-sided area), as in
%                      find_edfm_frac_rock.
%   PEDFMFracRock.Area is one-sided projected area.
%   PEDFMFracRock.AlphaID lists the faces whose alpha is set to 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EDFMFracRock, ...
          PEDFMFracRock ] = ...
          FracRockIntersectionMC(FracMesh, ...
                                 RockMesh, ...
                                 FracGeometry, ...
                                 Settings)

Settings = monte_carlo_settings(Settings);
PrevRng = rng;
RestoreRng = onCleanup(@() rng(PrevRng));
rng(Settings.seed, "twister");

SampledProjection = Settings.pedfm && isfield(Settings, "pedfm_projection") ...
                    && string(Settings.pedfm_projection) == "sampled";
CoreSettings = Settings;
CoreSettings.pedfm = SampledProjection;
[edfm, pedfm] = find_pedfm_edfm_int_area_MC(FracMesh, RockMesh, ...
                                            FracGeometry, CoreSettings);

toDouble = @(c) cellfun(@double, c, 'UniformOutput', false);

EDFMFracRock = struct();
EDFMFracRock.RockID = toDouble(edfm.vol_ind);
EDFMFracRock.FracID = toDouble(edfm.frac_ind);
EDFMFracRock.Area = cellfun(@(a) 2*double(a), edfm.area, ...
                            'UniformOutput', false);

PEDFMFracRock = struct();
if SampledProjection
    PEDFMFracRock.RockID = toDouble(pedfm.vol_ind);
    PEDFMFracRock.FracID = toDouble(pedfm.frac_ind);
    PEDFMFracRock.Area = toDouble(pedfm.area);
    PEDFMFracRock.AlphaID = unique(double(pedfm.alpha_ind(:)));
elseif Settings.pedfm
    PEDFMFracRock = connector_projection(EDFMFracRock, RockMesh, ...
                                         FracGeometry);
end

end

function PEDFMFracRock = connector_projection(EDFMFracRock, RockMesh, ...
                                              FracGeometry)
nFrac = numel(EDFMFracRock.RockID);
RockID = cell(nFrac,1);
FracID = cell(nFrac,1);
Area = cell(nFrac,1);
AlphaID = cell(nFrac,1);
for iFrac = 1:nFrac
    RockIntersections = struct();
    RockIntersections.RockCellID = EDFMFracRock.RockID{iFrac}(:);
    RockIntersections.FracCellID = EDFMFracRock.FracID{iFrac}(:);
    RockIntersections.AreaEDFM = EDFMFracRock.Area{iFrac}(:);
    RockIntersections.Coplanar = double(coplanar_rock_cells( ...
        RockIntersections.RockCellID, RockMesh, FracGeometry(iFrac).eq));
    if isempty(RockIntersections.RockCellID)
        Connections = zeros(0,2);
        AreaPEDFM = zeros(0,1);
        AlphaFaceID = zeros(0,1);
    else
        [Connections, AreaPEDFM, AlphaFaceID] = find_pdefm_frac_rock( ...
            RockIntersections, RockMesh, FracGeometry(iFrac));
    end
    if isempty(Connections)
        Connections = zeros(0,2);
    end
    RockID{iFrac} = double(Connections(:,1));
    FracID{iFrac} = double(Connections(:,2));
    Area{iFrac} = double(AreaPEDFM(:));
    AlphaID{iFrac} = double(AlphaFaceID(:));
end
PEDFMFracRock.RockID = RockID;
PEDFMFracRock.FracID = FracID;
PEDFMFracRock.Area = Area;
PEDFMFracRock.AlphaID = unique(cell2mat(AlphaID));

end
