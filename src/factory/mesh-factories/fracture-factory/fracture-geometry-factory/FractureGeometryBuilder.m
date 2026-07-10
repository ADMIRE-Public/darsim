%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FractureGeometryBuilder
    %FractureBuilder created for quadrilateral fractures
    properties
        num_frac
        fracs
        info
    end

    methods
        function obj = FractureGeometryBuilder(blueprint)
            obj.num_frac = ...
                numel(cell2mat({blueprint.FracturesProperties.isActive})');
            obj.info = obj.build_info(blueprint);
            frac_tag = obj.classify_fracuture_building_method(blueprint);
            obj.fracs = obj.build_fractures(frac_tag, blueprint);
            obj.fracs = obj.add_frac_equation(obj.fracs);
        end
    end

    methods (Access = private)
        function [fracs] = add_frac_equation(obj,fracs)
            for index = 1:obj.num_frac
                fracs(index).eq = create_plane_eq(obj.fracs(index).center,...
                    obj.fracs(index).normal);
            end
        end
        function frac = build_fractures(obj, tag, blueprint)
            frac(obj.num_frac) = struct();
            build_frac_coords = ...
                {@obj.build_fracuture_coords_from_center_coord ...
                @obj.build_fracuture_coords_from_corner_coord};
            ii = 1;
            for index = ...
                    find(cell2mat({blueprint.FracturesProperties.isActive}))
                frac_builder = build_frac_coords{tag(index)+1};
                frac(ii).tag_reader = index;
                [frac(ii).coords, frac(ii).center] = ...
                    frac_builder(blueprint.FracturesProperties(index));
                [frac(ii).normal, frac(ii).area] = ...
                    frac_normal(frac(ii).coords);
                % the area of the fracture is counted twice: up/down
                frac(ii).area = 2*frac(ii).area;
                frac(ii).area = frac(ii).area;
                [frac(ii).vec_quad_edges, frac(ii).edges] = ...
                    create_boundaries(frac(ii).coords);
                frac(ii).Aperture = ...
                    blueprint.FracturesProperties(index).Aperture;
                frac(ii).Porosity = ...
                    blueprint.FracturesProperties(index).Porosity;
                frac(ii).Permeability = ...
                    blueprint.FracturesProperties(index).Permeability;
                frac(ii).ADM = blueprint.FracturesProperties(index).ADM;
                frac(ii).GridResRatio = ...
                    blueprint.FracturesProperties(index).GridResRatio;
                frac(ii).GridNumAlongL = ...
                    blueprint.FracturesProperties(index).GridNumAlongL;
                frac(ii).GridNumAlongW = ...
                    blueprint.FracturesProperties(index).GridNumAlongW;
                ii = ii + 1;
            end
        end
        function info = build_info(obj,blueprint)
            info = struct();
            info.Type = blueprint.Type;
            info.NumOfDimensions = blueprint.NumOfDimensions;
            info.Discretization = blueprint.Discretization;
            info.ADM = blueprint.ADM;
            info.tol = blueprint.compareAccuracy;
            info.Size = blueprint.Size;
            info.N = blueprint.N;
        end
        function tag = classify_fracuture_building_method(obj, blueprint)
            % 0 - Create Square Fracture from from Center Coord.
            % 1 - Create Square Fracture from Corner Coords.
            tag = zeros(obj.num_frac,1);
            for index = 1:size(blueprint.FracturesProperties,1)
                flag = ...
                    ~all(isnan(blueprint.FracturesProperties...
                    (index).CornerCoords(:)));
                tag(index) = flag;
            end
        end

        function [points,frac_center] = ...
                build_fracuture_coords_from_center_coord(obj, frac_blueprint)
            %Create 2D square fracture from from Center Coord.
            frac_center = frac_blueprint.CenterCoord;
            Angles = [frac_blueprint.RotAngleAlongZ ...
                frac_blueprint.RotAngleAlongL ...
                frac_blueprint.RotAngleAlongW];
            points = create_square_frac(frac_center, ...
                frac_blueprint.Length, frac_blueprint.Width);
            points = rotate_points(frac_center, Angles, points);

        end

        function [coords, frac_center] = ...
                build_fracuture_coords_from_corner_coord(obj, frac_blueprint)
            % Create 2D square fracture from Corner Coords.
            % Assumes the CornerCoords is given in order
            points3 = frac_blueprint.CornerCoords;
            A = points3(1,:);
            B = points3(2,:);
            C = points3(3,:);
            AB =  A - B;
            CB =  C - B;
            AD = AB+ CB;
            D = AD + B;
            coords =[A;B;C;D];
            frac_center = mean(coords,1);
        end

    end
end

function points = create_square_frac(PointM, length, width)
half_length = 0.5*length;
half_width = 0.5*width;
PointA = PointM + [-half_length, half_width, 0];
PointB = PointM + [half_length, half_width, 0];
PointC = PointM + [half_length, -half_width, 0];
PointD = PointM + [-half_length, -half_width, 0];
points = [PointA;PointB;PointC;PointD];
end

function rot_points = rotate_points(center,angles, varargin)
points = cell2mat(reshape(varargin,[],1));
npoints = size(points,1);

A = points(1,:);
B = points(2,:);
C = points(3,:);
D = points(4,:);

RotZ_deg = angles(1);
RotL_deg = angles(2);
RotW_deg = angles(3);

Zaxis = [0 0 1];
A = rot_around_axis(A, center, Zaxis,  RotZ_deg);
B = rot_around_axis(B, center, Zaxis,  RotZ_deg);
C = rot_around_axis(C, center, Zaxis,  RotZ_deg);
D = rot_around_axis(D, center, Zaxis,  RotZ_deg);

BMC = 0.5*(B + C);
DMA = 0.5*(D + A);
Laxis = DMA - BMC;
A = rot_around_axis(A, center, Laxis, RotL_deg);
B = rot_around_axis(B, center, Laxis, RotL_deg);
C = rot_around_axis(C, center, Laxis, RotL_deg);
D = rot_around_axis(D, center, Laxis, RotL_deg);

AMB = 0.5*(A + B);
CMD = 0.5*(C + D);
Waxis = CMD - AMB;  
A = rot_around_axis(A, center, Waxis, RotW_deg);
B = rot_around_axis(B, center, Waxis, RotW_deg);
C = rot_around_axis(C, center, Waxis, RotW_deg);
D = rot_around_axis(D, center, Waxis, RotW_deg);

rot_points = [A; B; C; D];
end

function P = rot_around_axis(P, O, axis, ang_deg)
    v = axis(:); nv = norm(v);
    if nv < 1e-15 || abs(ang_deg) < 1e-15, return; end
    n = (v/nv).';                 % 1x3
    theta = deg2rad(ang_deg);
    c = cos(theta); s = sin(theta);
    nx=n(1); ny=n(2); nz=n(3);
    R = [nx*nx+(1-nx*nx)*c, nx*ny*(1-c)-nz*s, nx*nz*(1-c)+ny*s; ...
         nx*ny*(1-c)+nz*s, ny*ny+(1-ny*ny)*c, ny*nz*(1-c)-nx*s; ...
         nx*nz*(1-c)-ny*s, ny*nz*(1-c)+nx*s, nz*nz+(1-nz*nz)*c];
    P = (P - O) * R.' + O;
end


function [normal, area] = frac_normal(coords)
vec1 = coords(2,:) - coords(1,:);
vec2 = coords(3,:) - coords(1,:);
normal = cross(vec1,vec2);
area =  norm(normal);
normal = area .\ normal;
end

function [vec_edges, edge_length ] = create_boundaries(coords)
vec_edges = circshift(coords,-1) - coords;
edge_length = vecnorm(vec_edges,2,2);
end

function eq = create_plane_eq(coords,normal)
% Ax+By+Cz+D = 0
eq = struct();
eq.A = normal(1);
eq.B = normal(2);
eq.C = normal(3);
eq.D = -sum(normal.*coords(1,:));
end
