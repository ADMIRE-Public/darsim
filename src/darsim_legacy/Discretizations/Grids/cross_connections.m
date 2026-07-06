classdef cross_connections < handle
    % Contains the information of non-neighboring connectivities between different media (fracture0matrix and fracture-fracture)
    properties
        Cells
        CI
        T_Geo
        T_Geo_Cond
        UpWind
        U_Geo
    end
    methods
        function obj = cross_connections()
            
        end
        function CrossConnectionsUpWind(obj, P, Depth, rho, g, N, I)
            [~, n_phases] = size(obj.UpWind);
            Depth = Depth * 0;
            depthDiff    = (Depth(obj.Cells) - Depth(I)) ;
            for i = 1 : n_phases
                rhoint = (rho(obj.Cells,i) + rho(I,i) ) / 2 ;
                obj.UpWind(:,i) = ((P(obj.Cells) - P(I) -  rhoint .* g .* depthDiff) >= 0);
                obj.U_Geo(:,i) = obj.T_Geo .* (P(obj.Cells) - P(I) - rhoint .* g .* depthDiff);
            end
        end
    end
    
end

