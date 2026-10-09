%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% Calculating the physical properties of EDFM and pEDFM connections between
% fracture-rock cells and fracture-fracture cells.
%
% Input:
%   rockmesh:       struct
%                   Contains relevant rock mesh information
%   fracmesh:       struct
%                   Contains relevant fracture mesh information
%   geometry:       struct
%                   Contains relevant fracture geometry and properties
%   fracmeshint:    cell
%                   Contains fracture-fracture mesh intersection information
%   edfm_geo:       struct
%                   Contains geometric information for EDFM connections
%   pedfm_geo:      struct
%                   Contains geometric information for pEDFM connections
%   avg_d_rock:     cell
%                   Contains average distances between fracture and rock
%                   cells
%   avg_d_frac:     cell
%                   Contains average distances between fracture cells and
%                   fracture-fracture intersection lines
%
% Output:
%   pedfm:          struct
%                   Contains physical properties for intra-fracture,
%                   fracture-fracture, EDFM, and pEDFM connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pedfm] = ...
         compute_physical_prop_pedfm(rockmesh, ...
                                     fracmesh, ...
                                     geometry, ...
                                     fracmeshint, ...
                                     edfm_geo, ...
                                     pedfm_geo, ...
                                     avg_d_rock, ...
                                     avg_d_frac)

% Initialize output structure
pedfm = struct();

% Number of fractures
nFrac = size(fracmesh,1);

% Calculate fracture-fracture EDFM connection properties
FracFracTrans = cell(nFrac,nFrac);

for ii = 1:size(fracmeshint,1)

    for jj = ii+1:size(fracmeshint,2)

        if ~isempty(fracmeshint{ii,jj})

            % Initialize physical properties
            Physical = struct();

            % Fracture cells connected through fracture intersection
            LocalFrac = fracmeshint{ii,jj};
            LocalInd1 = LocalFrac(:,1);
            LocalInd2 = LocalFrac(:,2);

            % Length of fracture-fracture intersection within each
            % connected pair of fracture cells
            IntDist = LocalFrac(:,3);

            % Average distance from fracture cells to intersection line
            AuxCell = avg_d_frac{ii,jj};
            D1 = AuxCell{1};
            D2 = AuxCell{2};

            % Fracture permeabilities
            K1 = geometry(ii).Permeability;
            K2 = geometry(jj).Permeability;

            % Harmonic average permeability between fractures
            Physical.K = (2*K1*K2)./(K1+K2);

            % Map average distances to each fracture-fracture connection
            LocalD1 = apply_filter(LocalInd1, D1, unique(LocalInd1));
            LocalD2 = apply_filter(LocalInd2, D2, unique(LocalInd2));

            Physical.ind1 = LocalInd1;
            Physical.ind2 = LocalInd2;

            % Connectivity index from each fracture to intersection line
            CI1 = 2*IntDist./LocalD1;
            CI2 = 2*IntDist./LocalD2;

            % Harmonic average fracture-fracture connectivity index
            CIHarm = 2*(CI1.*CI2)./(CI1+CI2);
            Physical.CI = CIHarm;

            FracFracTrans{ii,jj} = Physical;

        end

    end

end

pedfm.fracfrac = FracFracTrans;

% Calculate fracture-rock EDFM connection properties
FracRockEDFMTrans = cell(nFrac,1);

for iFrac = 1:nFrac

    % Initialize physical properties
    Physical = struct();

    % Rock and fracture cell indices of EDFM connections
    VolInd = edfm_geo.vol_ind{iFrac};
    FracInd = edfm_geo.frac_ind{iFrac};

    % Fracture-rock intersection area
    SharedFracArea = edfm_geo.area{iFrac};

    % Average distance between fracture and rock cells
    LocalDInd = avg_d_rock{iFrac,1};
    LocalDValue = avg_d_rock{iFrac,2};
    DVolInd = apply_filter(VolInd, LocalDValue, LocalDInd);

    % EDFM connectivity index
    CI = SharedFracArea./DVolInd;

    % Rock and fracture permeability
    K1 = rockmesh.props.perm(VolInd);
    K2 = geometry(iFrac).Permeability;
    K2 = repmat(K2,size(K1,1),1);

    % Characteristic rock cell length and fracture aperture
    L = rockmesh.element.volume(VolInd).^(1/3);
    W = geometry(iFrac).Aperture;

    % Store EDFM connection properties
    Physical.vol_ind = VolInd;
    Physical.frac_ind = FracInd;

    % Effective permeability between rock and fracture cells
    Physical.K = (L+W).*(K1.*K2)./(W.*K1+L.*K2);
    Physical.CI = CI;

    FracRockEDFMTrans{iFrac} = Physical;

end

pedfm.fracrock_edfm = FracRockEDFMTrans;

% Initialize pEDFM properties
pedfm.pedfm_on = false;
pedfm.alpha = zeros(size(rockmesh.face.inner.area));

% Calculate fracture-rock pEDFM connection properties
if ~isempty(fieldnames(pedfm_geo))

    pedfm.pedfm_on = true;
    FracRockPEDFMTrans = cell(size(fracmeshint,1),1);

    for iFrac = 1:nFrac

        % Initialize physical properties
        Physical = struct();

        % Rock and fracture cell indices of pEDFM connections
        VolInd = pedfm_geo.vol_ind{iFrac};
        FracInd = pedfm_geo.frac_ind{iFrac};

        % Projected fracture-rock intersection area
        SharedFracArea = pedfm_geo.area{iFrac};

        % Average distance between fracture and rock cells
        LocalDInd = avg_d_rock{iFrac,1};
        LocalDValue = avg_d_rock{iFrac,2};
        DVolInd = ...
            apply_pdefm_filter(VolInd, LocalDInd, LocalDValue);

        % pEDFM connectivity index
        CI = SharedFracArea./DVolInd;

        % Rock and fracture permeability
        K1 = rockmesh.props.perm(VolInd);
        K2 = geometry(iFrac).Permeability;
        K2 = repmat(K2,size(K1,1),1);

        % Characteristic rock cell length and fracture aperture
        L = rockmesh.element.volume(VolInd).^(1/3);
        W = geometry(iFrac).Aperture;

        % Effective permeability between rock and fracture cells
        Physical.K = (L+W).*(K1.*K2)./(W.*K1+L.*K2);
        Physical.CI = CI;
        Physical.vol_ind = VolInd;
        Physical.frac_ind = FracInd;

        FracRockPEDFMTrans{iFrac} = Physical;

    end

    pedfm.fracrock_pedfm = FracRockPEDFMTrans;

end

% Calculate pEDFM alpha values
if ~isStructEmpty(pedfm_geo)

    % Fraction of face area required to block matrix-matrix connection
    BlockTol = 0.1;

    % Areas of all rock mesh faces
    AreaTotal = [rockmesh.face.inner.area; ...
                 rockmesh.face.bound.area];

    % Number of rock mesh faces
    nFace = numel(AreaTotal);

    % Sum projected fracture area associated with each rock mesh face
    FaceAreaLoss = accumarray(pedfm_geo.alpha_ind(:), ...
                              pedfm_geo.alpha_area(:), ...
                              [nFace 1], @sum, 0);

    % Effective face area remaining after fracture projection
    AreaEffective = AreaTotal - FaceAreaLoss;
    AreaEffective(AreaEffective < 0) = 0;

    % Completely block faces exceeding the specified area threshold
    BlockedFrac = FaceAreaLoss./AreaTotal;
    AreaEffective(BlockedFrac > BlockTol) = 0;

    % Fraction of rock mesh face occupied by fracture projection
    AlphaGeom = FaceAreaLoss./AreaTotal;
    AlphaGeom(AlphaGeom > 1) = 1;

    % Store alpha values
    pedfm.alpha = AlphaGeom;
    pedfm.alpha(AreaTotal == 0) = 0;

end

end

function [filter] = apply_filter(ind, value, allInd)

    % Create mapping from cell indices to associated values
    allInd = double(allInd);
    TransfMatrix = sparse(allInd, ones(length(allInd),1), value);

    % Extract values corresponding to requested cell indices
    filter = full(TransfMatrix(ind));

end

function [filter] = apply_pdefm_filter(input, ind, value)

    % Create mapping from rock cell indices to associated distance values
    ind = double(ind);
    OneInd = ones(size(ind));
    MaxValue = max(ind);

    TransfMatrix = sparse(ind, OneInd, value, MaxValue, 1);

    % Extract values corresponding to pEDFM rock cell indices
    filter = full(TransfMatrix(input));

end