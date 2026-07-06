function [FracturesGrid, CrossConnections] = legacy_scanfracture_data(obj,ReservoirGrid)
%warning("fracture data read from txt file is transformed into matlab data here")
fprintf('This simulation uses input from "%s" fracture generator.\n', obj.SimulationInput.FracturesProperties.FractureModelType );
if ~strcmp( obj.SimulationInput.FracturesProperties.DiscretizationType, obj.SimulationInput.ReservoirProperties.Discretization )
    error('The discretization type used in "Fracture_Output" file (%s) does not match the one in the simulation input file (%s)!', ...
        obj.SimulationInput.FracturesProperties.DiscretizationType, obj.SimulationInput.ReservoirProperties.Discretization );
end
NumOfFrac = obj.SimulationInput.FracturesProperties.NumOfFrac;
fprintf('Extracting data from %02d fractures ...\n', NumOfFrac);
FracturesGrid = fractures_grid(NumOfFrac);
if strcmp(obj.SimulationInput.ReservoirProperties.Discretization,'TetrahedralGrid')
    if ( ReservoirGrid.N ~= obj.SimulationInput.FracturesProperties.Reservoir.N )
        error('The number of grid cells for reservoir in the "Fracture_Output" file does not match the simulation input file!');
    end
else
    if ( obj.SimulationInput.ReservoirProperties.size(1) ~= obj.SimulationInput.FracturesProperties.Reservoir.Lx ) || ...
            ( obj.SimulationInput.ReservoirProperties.size(2) ~= obj.SimulationInput.FracturesProperties.Reservoir.Ly ) || ...
            ( obj.SimulationInput.ReservoirProperties.size(3) ~= obj.SimulationInput.FracturesProperties.Reservoir.Lz )
        warning('DARSim Warning: The dimension of reservoir in the "Fracture_Output" file does not match the simulation input file. This may cause inconsistency!');
    end

    if ( ReservoirGrid.Nx ~= obj.SimulationInput.FracturesProperties.Reservoir.Nx ) || ...
            ( ReservoirGrid.Ny ~= obj.SimulationInput.FracturesProperties.Reservoir.Ny ) || ...
            ( ReservoirGrid.Nz ~= obj.SimulationInput.FracturesProperties.Reservoir.Nz )
        error('The number of grid cells for reservoir in the "Fracture_Output" file does not match the simulation input file!');
    end

end

fileID = fopen(obj.SimulationInput.FracturesProperties.FilePath, 'r');
fprintf('Reading fracture file ...');
matrix = textscan(fileID, '%s', 'Delimiter', '\n');
FractureMatrix = matrix{1};
fclose(fileID);
temp = strfind(FractureMatrix, 'PROPERTIES');
frac_index = find(~cellfun('isempty', temp));
Nf = zeros(NumOfFrac,1);
%warning("creation of cartesian grid in each fracture is done here")
for f = 1 : NumOfFrac
    % Creat cartesian grid in each fracture
    frac_info_split = strsplit( FractureMatrix{frac_index(f)} , ' ' );
    grid_temp = strsplit(frac_info_split{8}, 'x');
    nx = str2double(grid_temp{1});
    ny = str2double(grid_temp{2});
    nz = 1;
    FractureGrid = cartesian_grid([nx;ny;nz]);

    % Add fracture grid coordinates (it's for plotting purposes)
    temp = strfind(FractureMatrix, 'GRID_COORDS_X');  frac_grid_coords_x = find(~cellfun('isempty', temp));
    temp = strfind(FractureMatrix, 'GRID_COORDS_Y');  frac_grid_coords_y = find(~cellfun('isempty', temp));
    temp = strfind(FractureMatrix, 'GRID_COORDS_Z');  frac_grid_coords_z = find(~cellfun('isempty', temp));

    frac_grid_coords_x_split = str2double( strsplit(FractureMatrix{frac_grid_coords_x(f)},' ') );  frac_grid_coords_x_split(1) = [];
    frac_grid_coords_y_split = str2double( strsplit(FractureMatrix{frac_grid_coords_y(f)},' ') );  frac_grid_coords_y_split(1) = [];
    frac_grid_coords_z_split = str2double( strsplit(FractureMatrix{frac_grid_coords_z(f)},' ') );  frac_grid_coords_z_split(1) = [];
    FractureGrid.GridCoords = [ frac_grid_coords_x_split' , frac_grid_coords_y_split' , frac_grid_coords_z_split' ];
    FracturesGrid.AddGrid(FractureGrid, f);
    Nf(f) = nx*ny*nz;
end

% Reading the non-neighboring connectivities of frac-frac and frac-matrix
NumOfAllFracGrids = obj.SimulationInput.FracturesProperties.NumOfAllFracGrids;
CrossConnections = cross_connections();

temp = strfind(FractureMatrix, '/');
endOfSection = find(~cellfun('isempty', temp));

temp = strfind(FractureMatrix, 'FRACCELL');
fracCellIndeces = find(~cellfun('isempty', temp));
fracCellIndeces = [fracCellIndeces; endOfSection(end)];

temp = strfind(FractureMatrix, 'ROCK_CONN_EDFM');
rockConnIndeces_EDFM = find(~cellfun('isempty', temp));

temp = strfind(FractureMatrix, 'ROCK_CONN_pEDFM');
rockConnIndeces_pEDFM = find(~cellfun('isempty', temp));

temp = strfind(FractureMatrix, 'FRAC_CONN_EDFM');
fracConnIndeces_EDFM = find(~cellfun('isempty', temp));

temp = strfind(FractureMatrix, 'FRAC_CONN_pEDFM');
fracConnIndeces_pEDFM = find(~cellfun('isempty', temp));

% Number of phases (useful to define size of some objects)
n_phases = obj.SimulationInput.FluidProperties.NofPhases;
fprintf('---> Fracture ');

Nm = ReservoirGrid.N;
for f = 1 : NumOfFrac
    if (f>1),  fprintf(repmat('\b', 1, 9+27));  end
    fprintf('%04d/%04d',f,NumOfFrac);
    % looping over all global fracture cells
    fprintf(' ---> Grid cell ');
    for If = 1:Nf(f)
        if (If>1),  fprintf(repmat('\b', 1, 11));  end
        fprintf('%05d/%05d',If,Nf(f));
        IfGlobal = sum(Nf(1:f-1))+If;
        fracCellInfo_split = strsplit( FractureMatrix{fracCellIndeces(IfGlobal)} , {' ','	'} );
        fracCellInd = str2double(fracCellInfo_split{2})+1;
        if If~= fracCellInd,  error('The Indexing for fracture %02d cell %05d is not correct!',f,If);  end
        NumOfRockConn_EDFM  = str2double(fracCellInfo_split{3});
        NumOfRockConn_pEDFM = str2double(fracCellInfo_split{4});
        NumOfFracConn_EDFM  = str2double(fracCellInfo_split{5});
        NumOfFracConn_pEDFM = str2double(fracCellInfo_split{6});

        Counter=0;

        % Adding fracture-matrix EDFM connectivities
        CrossConnections(IfGlobal,1).Cells = [];
        rockConnInd_EDFM = rockConnIndeces_EDFM( rockConnIndeces_EDFM > fracCellIndeces(IfGlobal) );
        rockConnInd_EDFM = rockConnInd_EDFM( rockConnInd_EDFM < fracCellIndeces(IfGlobal+1) );
        if length(rockConnInd_EDFM) ~= NumOfRockConn_EDFM
            error('The number of rock EDFM connectivities for fracture %02d cell %05d is not correct!',f,If);
        end
        for m = 1 : NumOfRockConn_EDFM
            rockConnInfo = strsplit( FractureMatrix{rockConnInd_EDFM(m)} ,  {' ','	'} );
            Im = str2double( rockConnInfo{2} ) + 1;
            CI = str2double( rockConnInfo{3} );
            Counter = Counter + 1;
            CrossConnections( IfGlobal , 1 ).Cells( Counter ,1 ) = Im;
            CrossConnections( IfGlobal , 1 ).CI   ( Counter ,1 ) = CI;
        end

        % Adding fracture-matrix pEDFM connectivities
        if strcmp(obj.SimulationInput.FracturesProperties.FractureModelType,'pEDFM')
            rockConnInd_pEDFM = rockConnIndeces_pEDFM( rockConnIndeces_pEDFM > fracCellIndeces(IfGlobal) );
            rockConnInd_pEDFM = rockConnInd_pEDFM( rockConnInd_pEDFM < fracCellIndeces(IfGlobal+1) );
            if length(rockConnInd_pEDFM) ~= NumOfRockConn_pEDFM
                error('The number of rock pEDFM connectivities for fracture %02d cell %05d is not correct!',f,If);
            end
            for m = 1 : NumOfRockConn_pEDFM
                rockConnInfo = strsplit( FractureMatrix{rockConnInd_pEDFM(m)} ,  {' ','	'} );
                Im = str2double( rockConnInfo{2} ) + 1;
                CI = str2double( rockConnInfo{3} );
                if ismember( Im, CrossConnections( IfGlobal , 1 ).Cells )
                    ind = find(CrossConnections( IfGlobal , 1 ).Cells == Im);
                    CrossConnections( IfGlobal , 1 ).CI( ind ,1 ) = CrossConnections( IfGlobal , 1 ).CI( ind ,1 ) + CI;
                else
                    Counter = Counter + 1;
                    CrossConnections( IfGlobal , 1 ).Cells( Counter ,1 ) = Im;
                    CrossConnections( IfGlobal , 1 ).CI   ( Counter ,1 ) = CI;
                end
            end
        end

        % Adding fracture-fracture EDFM connectivities
        fracConnInd_EDFM = fracConnIndeces_EDFM( fracConnIndeces_EDFM > fracCellIndeces(IfGlobal) );
        fracConnInd_EDFM = fracConnInd_EDFM( fracConnInd_EDFM < fracCellIndeces(IfGlobal+1) );
        if length(fracConnInd_EDFM) ~= NumOfFracConn_EDFM
            error('The number of frac EDFM connectivities for fracture %02d cell %05d is not correct!',f,If);
        end
        for m = 1 : NumOfFracConn_EDFM
            fracConnInfo = strsplit( FractureMatrix{fracConnInd_EDFM(m)} ,  {' ','	'} );
            g = str2double( fracConnInfo{2} ) + 1;
            ig = str2double( fracConnInfo{3} ) + 1;
            Ig = sum( Nf(1:g-1) )+ ig;
            CI = str2double( fracConnInfo{4} );
            if Ig > IfGlobal
                Counter = Counter + 1;
                CrossConnections( IfGlobal , 1 ).Cells( Counter ,1 ) = Nm+Ig;
                CrossConnections( IfGlobal , 1 ).CI   ( Counter ,1 ) = CI;
            end
        end

        % Adding fracture-fracture pEDFM connectivities
        if strcmp(obj.SimulationInput.FracturesProperties.FractureModelType,'pEDFM')
            fracConnInd_pEDFM = fracConnIndeces_pEDFM( fracConnIndeces_pEDFM > fracCellIndeces(IfGlobal) );
            fracConnInd_pEDFM = fracConnInd_pEDFM( fracConnInd_pEDFM < fracCellIndeces(IfGlobal+1) );
            for m = 1 : NumOfFracConn_pEDFM
                fracConnInfo = strsplit( FractureMatrix{fracConnInd_pEDFM(m)} ,  {' ','	'} );
                g = str2double( fracConnInfo{2} ) + 1;
                ig = str2double( fracConnInfo{3} ) + 1;
                Ig = sum( Nf(1:g-1) )+ ig;
                CI = str2double( fracConnInfo{4} );
                if Ig > IfGlobal
                    if ismember( Nm+Ig, CrossConnections( IfGlobal , 1 ).Cells )
                        ind = find(CrossConnections( IfGlobal , 1 ).Cells == Nm+Ig);
                        CrossConnections( IfGlobal , 1 ).CI( ind ,1 ) = CrossConnections( IfGlobal , 1 ).CI( ind ,1 ) + CI;
                    else
                        Counter = Counter + 1;
                        CrossConnections( IfGlobal , 1 ).Cells( Counter ,1 ) = Nm+Ig;
                        CrossConnections( IfGlobal , 1 ).CI   ( Counter ,1 ) = CI;
                    end
                end
            end
        end

        % Initializing CrossConnections
        [~, sort_ind] = sort( CrossConnections(IfGlobal,1).Cells );
        CrossConnections(IfGlobal,1).Cells  = CrossConnections(IfGlobal,1).Cells(sort_ind);
        CrossConnections(IfGlobal,1).CI     = CrossConnections(IfGlobal,1).CI(sort_ind);
        CrossConnections(IfGlobal,1).UpWind = zeros(length(CrossConnections(IfGlobal,1).Cells), n_phases);
        CrossConnections(IfGlobal,1).U_Geo  = zeros(length(CrossConnections(IfGlobal,1).Cells), n_phases);
        CrossConnections(IfGlobal,1).T_Geo  = zeros(  size(CrossConnections(IfGlobal,1).CI)  );
        CrossConnections(IfGlobal,1).T_Geo_Cond  = zeros(  size(CrossConnections(IfGlobal,1).CI)  );
    end
end

% Reading the pEDFM alpha corrections for reservoir transmissiblities
if strcmp(obj.SimulationInput.FracturesProperties.FractureModelType,'pEDFM')

    % 1. Reservoir
    switch obj.SimulationInput.FracturesProperties.DiscretizationType
        case('CartesianGrid')
            pEDFM_alpha_Tx = zeros(ReservoirGrid.Nx+1,ReservoirGrid.Ny,ReservoirGrid.Nz);
            pEDFM_alpha_Ty = zeros(ReservoirGrid.Nx,ReservoirGrid.Ny+1,ReservoirGrid.Nz);
            pEDFM_alpha_Tz = zeros(ReservoirGrid.Nx,ReservoirGrid.Ny,ReservoirGrid.Nz+1);
            % pEDFM_alpha_Tx
            temp = strfind(FractureMatrix, 'ROCK_ALPHA_TX');
            pEDFM_alpha_Tx_index = find(~cellfun('isempty', temp));
            for t = 1 : length(pEDFM_alpha_Tx_index)
                pEDFM_alpha_Tx_Split = strsplit(FractureMatrix{pEDFM_alpha_Tx_index(t)},' ');
                i     = str2double(pEDFM_alpha_Tx_Split{2})+1;
                j     = str2double(pEDFM_alpha_Tx_Split{3})+1;
                k     = str2double(pEDFM_alpha_Tx_Split{4})+1;
                alpha = str2double(pEDFM_alpha_Tx_Split{5});
                pEDFM_alpha_Tx(i,j,k) = alpha;
            end
            % pEDFM_alpha_Ty
            temp = strfind(FractureMatrix, 'ROCK_ALPHA_TY');
            pEDFM_alpha_Ty_index = find(~cellfun('isempty', temp));
            for t = 1 : length(pEDFM_alpha_Ty_index)
                pEDFM_alpha_Ty_Split = strsplit(FractureMatrix{pEDFM_alpha_Ty_index(t)},' ');
                i     = str2double(pEDFM_alpha_Ty_Split{2})+1;
                j     = str2double(pEDFM_alpha_Ty_Split{3})+1;
                k     = str2double(pEDFM_alpha_Ty_Split{4})+1;
                alpha = str2double(pEDFM_alpha_Ty_Split{5});
                pEDFM_alpha_Ty(i,j,k) = alpha;
            end
            % pEDFM_alpha_Tz
            temp = strfind(FractureMatrix, 'ROCK_ALPHA_TZ');
            pEDFM_alpha_Tz_index = find(~cellfun('isempty', temp));
            for t = 1 : length(pEDFM_alpha_Tz_index)
                pEDFM_alpha_Tz_Split = strsplit(FractureMatrix{pEDFM_alpha_Tz_index(t)},' ');
                i     = str2double(pEDFM_alpha_Tz_Split{2})+1;
                j     = str2double(pEDFM_alpha_Tz_Split{3})+1;
                k     = str2double(pEDFM_alpha_Tz_Split{4})+1;
                alpha = str2double(pEDFM_alpha_Tz_Split{5});
                pEDFM_alpha_Tz(i,j,k) = alpha;
            end
            ReservoirGrid.AddpEDFMCorrections( pEDFM_alpha_Tx , pEDFM_alpha_Ty , pEDFM_alpha_Tz );

        case('CornerPointGrid')
            pEDFM_alpha_T = zeros( length(ReservoirGrid.Trans) , 1 );
            temp = strfind(FractureMatrix, 'ROCK_pEDFM_ALPHA');
            pEDFM_alpha_T_index = find(~cellfun('isempty', temp));
            for t = 1 : length(pEDFM_alpha_T_index)
                pEDFM_alpha_T_Split = strsplit(FractureMatrix{pEDFM_alpha_T_index(t)},' ');
                i     = str2double(pEDFM_alpha_T_Split{2})+1;
                alpha = str2double(pEDFM_alpha_T_Split{3});
                pEDFM_alpha_T(i) = alpha;
            end
            ReservoirGrid.AddpEDFMCorrections( pEDFM_alpha_T );

        case('TetrahedralGrid')
            pEDFM_alpha_T = zeros( length(ReservoirGrid.Trans) , 1 );
            temp = strfind(FractureMatrix, 'ROCK_pEDFM_ALPHA');
            pEDFM_alpha_T_index = find(~cellfun('isempty', temp));
            for t = 1 : length(pEDFM_alpha_T_index)
                pEDFM_alpha_T_Split = strsplit(FractureMatrix{pEDFM_alpha_T_index(t)},' ');
                i     = str2double(pEDFM_alpha_T_Split{2})+1;
                alpha = str2double(pEDFM_alpha_T_Split{3});
                pEDFM_alpha_T(i) = alpha;
            end
            ReservoirGrid.AddpEDFMCorrections( pEDFM_alpha_T );

    end
    % 2. Fractures
    pEDFM_alpha_Tx = cell(FracturesGrid.Nfrac,1);
    pEDFM_alpha_Ty = cell(FracturesGrid.Nfrac,1);
    pEDFM_alpha_Tz = 0;
    for f = 1 : FracturesGrid.Nfrac
        pEDFM_alpha_Tx{f} = zeros(FracturesGrid.Grids(f).Nx+1,FracturesGrid.Grids(f).Ny);
        pEDFM_alpha_Ty{f} = zeros(FracturesGrid.Grids(f).Nx,FracturesGrid.Grids(f).Ny+1);
    end
    % Tx_alpha
    temp = strfind(FractureMatrix, 'FRAC_ALPHA_TX');
    pEDFM_alpha_Tx_index = find(~cellfun('isempty', temp));
    for t = 1 : length(pEDFM_alpha_Tx_index)
        pEDFM_alpha_Tx_Split = strsplit(FractureMatrix{pEDFM_alpha_Tx_index(t)},' ');
        f   = str2double(pEDFM_alpha_Tx_Split{2})+1;
        i_f = str2double(pEDFM_alpha_Tx_Split{3})+1;
        j_f = str2double(pEDFM_alpha_Tx_Split{4})+1;
        alpha = str2double(pEDFM_alpha_Tx_Split{5});
        pEDFM_alpha_Tx{f}(i_f,j_f) = alpha;
    end
    % Ty_alpha
    temp = strfind(FractureMatrix, 'FRAC_ALPHA_TY');
    pEDFM_alpha_Ty_index = find(~cellfun('isempty', temp));
    for t = 1 : length(pEDFM_alpha_Ty_index)
        pEDFM_alpha_Ty_Split = strsplit(FractureMatrix{pEDFM_alpha_Ty_index(t)},' ');
        f   = str2double(pEDFM_alpha_Ty_Split{2})+1;
        i_f = str2double(pEDFM_alpha_Ty_Split{3})+1;
        j_f = str2double(pEDFM_alpha_Ty_Split{4})+1;
        alpha = str2double(pEDFM_alpha_Ty_Split{5});
        pEDFM_alpha_Ty{f}(i_f,j_f) = alpha;
    end
    for f = 1 : FracturesGrid.Nfrac
        FracturesGrid.Grids(f).AddpEDFMCorrections(pEDFM_alpha_Tx{f},pEDFM_alpha_Ty{f},pEDFM_alpha_Tz)
    end
end

fprintf(' ---> Complete!\n');
end
