%  Tetrahedral grid
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Túlio de M. Cavalcante
% UFPE
%--------------------------------------------------------------------------
classdef TetrahedralGrid_Discretization_model < FS_Discretization_model
    properties
        TetrahedralGridData
    end
    methods
        function AddHarmonicPermeabilities(obj, Reservoir, Fractures)
            Nm = obj.ReservoirGrid.N;
            for If1Local = 1 : length(obj.CrossConnections)
                if If1Local == 21
                    warning("20")
                end

                %----------------------------------------------------------
                coords = obj.FracturesGrid.Grids.GridCoords;
                v1 = coords(2,:) - coords(1,:); v2 = coords(end,:) - coords(1,:);
                n = cross(v1,v2)/norm(cross(v1,v2));
                %----------------------------------------------------------
                thetax = Reservoir.KRotAngle(1); thetay = Reservoir.KRotAngle(2); thetaz = Reservoir.KRotAngle(3);
                Rx = [1 0 0; 0 cos(thetax) -sin(thetax); 0 sin(thetax) cos(thetax)];
                Ry = [cos(thetay) 0 sin(thetay); 0 1 0; -sin(thetay) 0 cos(thetay)];
                Rz = [cos(thetaz) -sin(thetaz) 0; sin(thetaz) cos(thetaz) 0; 0 0 1];
                %----------------------------------------------------------       
                Cells = obj.CrossConnections(If1Local).Cells;
                CI = obj.CrossConnections(If1Local).CI;
                If1Global = Nm+If1Local; % Global index of this fracture cell;
                Ind_frac1_Local = obj.Index_Global_to_Local(If1Global);
                f1 = Ind_frac1_Local.f;
                g1 = Ind_frac1_Local.g;
                indices_m = Cells( Cells <= Nm );
                %----------------------------------------------------------     
                ResK = zeros(size(indices_m));
                for ik = 1:size(indices_m,1)
                    ResK(ik) = n*Rz'*Ry'*Rx'*diag(Reservoir.K(indices_m(ik),:))*Rx*Ry*Rz*n';
                end
                %----------------------------------------------------------                    
                if ~isempty(indices_m)
                    obj.CrossConnections(If1Local).T_Geo(1:length(indices_m)) = ...
                        CI(1:length(indices_m)) .* ( obj.ReservoirGrid.Volume(indices_m).^(1/3) + Fractures(f1).Thickness ) ./...
                        ( ( obj.ReservoirGrid.Volume(indices_m).^(1/3) ./ ResK ) + ( Fractures(f1).Thickness ./ Fractures(f1).K(g1,1) ) );
                end
                
                indices_f = Cells( Cells > Nm );
                if ~isempty(indices_f)
                    for n = 1:length(indices_f)
                        If2Global = indices_f(n);
                        If2Local = If2Global - Nm;
                        Ind_frac2_Local = obj.Index_Global_to_Local(indices_f(n));
                        f2 = Ind_frac2_Local.f;
                        g2 = Ind_frac2_Local.g;
                        obj.CrossConnections(If1Local).T_Geo(length(indices_m)+n) = ...
                            CI(length(indices_m)+n) * ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 + Fractures(f2).Thickness ) ./...
                              ( ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 ./ Fractures(f1).K(g1,1) ) + ( Fractures(f2).Thickness ./ Fractures(f2).K(g2,1) ) );
                    end
                end
            end
        end
        function AddHarmonicConductivities(obj, Reservoir, Fractures)
            Nm = obj.ReservoirGrid.N;
            for If1Local = 1 : length(obj.CrossConnections)
                Cells = obj.CrossConnections(If1Local).Cells;
                CI = obj.CrossConnections(If1Local).CI;
                If1Global = Nm+If1Local; % Global index of this fracture cell;
                Ind_frac1_Local = obj.Index_Global_to_Local(If1Global);
                f1 = Ind_frac1_Local.f;
                g1 = Ind_frac1_Local.g;
                
                indices_m = Cells( Cells <= Nm );
                if ~isempty(indices_m)
                    obj.CrossConnections(If1Local).T_Geo_Cond(1:length(indices_m)) = ...
                        CI(1:length(indices_m)) .* ( obj.ReservoirGrid.Volume(indices_m).^(1/3) + Fractures(f1).Thickness ) ./...
                        ( ( obj.ReservoirGrid.Volume(indices_m).^(1/3) ./ Reservoir.K_Cond_eff(indices_m,1) ) + ( Fractures(f1).Thickness ./ Fractures(f1).K_Cond_eff(g1,1) ) );
                end
                
                indices_f = Cells( Cells > Nm );
                if ~isempty(indices_f)
                    for n = 1:length(indices_f)
                        If2Global = indices_f(n);
                        If2Local = If2Global - Nm;
                        Ind_frac2_Local = obj.Index_Global_to_Local(indices_f(n));
                        f2 = Ind_frac2_Local.f;
                        g2 = Ind_frac2_Local.g;
                        obj.CrossConnections(If1Local).T_Geo_Cond(length(indices_m)+n) = ...
                            CI(length(indices_m)+n) * ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 + Fractures(f2).Thickness ) ./...
                              ( ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 ./ Fractures(f1).K_Cond_eff(g1,1) ) + ( Fractures(f2).Thickness ./ Fractures(f2).K_Cond_eff(g2,1) ) );
                    end
                end
            end
        end
        function ObtainPerforatedCellsBasedOnXYZList(obj, Well, Well_Type, w)
            nvert = size(obj.TetrahedralGridData.Vertex.Coordinates,1);
            vert = obj.TetrahedralGridData.Vertex.Coordinates(:,1:2);
            influradius = mean(obj.TetrahedralGridData.Cell.Volume)^(1/3);
            wellpos = repmat(Well.Coordinate.Value(1,1:2),[nvert 1]);
            posvec = vert - wellpos; dist = sqrt(sum(posvec.*posvec,2));
            vertinwell = dist < influradius;
            [ ~, cellinwell ] = find(obj.TetrahedralGridData.Vertex.VertCellMatrix(find(vertinwell==1),:)==1);
            cellinwell = unique(cellinwell);
            wellcell = cellinwell( find( obj.TetrahedralGridData.Cell.Centroid(cellinwell,3) == max(obj.TetrahedralGridData.Cell.Centroid(cellinwell,3)) ) );
            cellperto = obj.TetrahedralGridData.Cell.Centroid(wellcell,1:2);
            wellpos = repmat(Well.Coordinate.Value(1,1:2),[size(cellperto,1) 1]);
            posvec = cellperto - wellpos; dist = sqrt(sum(posvec.*posvec,2));
            Well.Cells = [ Well.Cells ; wellcell( find( dist == min(dist)) ) ];
            Well.ResizeObjects(length(Well.Cells));
        end
        function ObtainPerforatedCellsBasedOnCellIndexList(obj, Well, Well_Type, w)
            Well.Cells = zeros(length(Well.Coordinate.Value),1);
            Well.Cells(:,1) = Well.Coordinate.Value;
            Well.Cells = unique(Well.Cells);
            Well.ResizeObjects(length(Well.Cells));
        end
    end
end