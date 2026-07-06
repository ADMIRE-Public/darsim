%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef UnstructuredGrid_Discretization_model < FS_Discretization_model
    properties
        UnstructuredGridData
    end
    methods
        function AddHarmonicPermeabilities(obj, Reservoir, Fractures)
            Nm = obj.ReservoirGrid.N;
            for If1Local = 1 : length(obj.CrossConnections)
                If1Global = Nm+If1Local; % Global index of this fracture cell;
                Ind_frac1_Local = obj.Index_Global_to_Local(If1Global);
                f1 = Ind_frac1_Local.f;
                g1 = Ind_frac1_Local.g;

                %----------------------------------------------------------
                %coords = obj.FracturesGrid.Grids.GridCoords(f1);
                % v1 = coords(2,:) - coords(1,:); v2 = coords(end,:) - coords(1,:);
                n = obj.FracturesGrid.dfrac.geometry(f1).normal;
                %----------------------------------------------------------
                thetax = Reservoir.KRotAngle(1); thetay = Reservoir.KRotAngle(2); thetaz = Reservoir.KRotAngle(3);
                Rx = [1 0 0; 0 cos(thetax) -sin(thetax); 0 sin(thetax) cos(thetax)];
                Ry = [cos(thetay) 0 sin(thetay); 0 1 0; -sin(thetay) 0 cos(thetay)];
                Rz = [cos(thetaz) -sin(thetaz) 0; sin(thetaz) cos(thetaz) 0; 0 0 1];
                %----------------------------------------------------------       
                Cells = double(obj.CrossConnections(If1Local).Cells);
                CI = obj.CrossConnections(If1Local).CI;
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
                        mesh1 = obj.FracturesGrid.Grids(f1).fracmesh;
                        all_d = [mesh1.faces.inner.length; mesh1.faces.bound.length];
                        faces_el1 = mesh1.elements.faces(g1,:);
                        w1 = Fractures(f1).Thickness;
                        K1 = Fractures(f1).K;
                        K1 = K1(g1);
                        L1 = mean(all_d(faces_el1));
                        %mesh2 = obj.FracturesGrid.Grids(f2).fracmesh;
                        w2 = Fractures(f2).Thickness;
                        K2 = Fractures(f2).K;
                        K2 = K2(g2);
                        Keq = (L1 + w1)/((L1/K1) + (w2/K2));
                        obj.CrossConnections(If1Local).T_Geo(length(indices_m)+n) = ...
                            CI(length(indices_m)+n) * Keq;
                        % obj.CrossConnections(If1Local).T_Geo(length(indices_m)+n) = ...
                        %     CI(length(indices_m)+n) * ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 + Fractures(f2).Thickness ) ./...
                        %       ( ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 ./ Fractures(f1).K(g1,1) ) + ( Fractures(f2).Thickness ./ Fractures(f2).K(g2,1) ) );
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
            P1 = Well.Coordinate.Value(1,:);
            P2 = Well.Coordinate.Value(2,:);
            radius = mean(obj.UnstructuredGridData.element.volume)^(1/3);
            test_points = obj.UnstructuredGridData.element.centroid;
            in_nodes = points_in_cylinder(P1',P2',radius,test_points')';
            if sum(in_nodes) == 0
                radius = radius*1.4;
                in_nodes = points_in_cylinder(P1',P2',radius,test_points')';
            end

            Well.Cells = find(in_nodes);
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
