%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% Monte Carlo alternative to FracRockIntersection (EDFM_METHOD
% monte_carlo). Test points sampled on every fracture cell are assigned to
% the rock cells that contain them; hit ratios give the intersection areas
% (EDFM) and the share of points projected inside neighbouring faces gives
% the projected areas (pEDFM).
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

[edfm, pedfm] = find_pedfm_edfm_int_area_MC(FracMesh, RockMesh, ...
                                            FracGeometry, Settings);

toDouble = @(c) cellfun(@double, c, 'UniformOutput', false);

EDFMFracRock = struct();
EDFMFracRock.RockID = toDouble(edfm.vol_ind);
EDFMFracRock.FracID = toDouble(edfm.frac_ind);
EDFMFracRock.Area = cellfun(@(a) 2*double(a), edfm.area, ...
                            'UniformOutput', false);

PEDFMFracRock = struct();
if Settings.pedfm
    PEDFMFracRock.RockID = toDouble(pedfm.vol_ind);
    PEDFMFracRock.FracID = toDouble(pedfm.frac_ind);
    PEDFMFracRock.Area = toDouble(pedfm.area);
    PEDFMFracRock.AlphaID = unique(double(pedfm.alpha_ind(:)));
end

end
