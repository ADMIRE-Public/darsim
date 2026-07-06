%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function node_classification = ...
                     classify_nodes_2d(p1,p2, coords)
%%
%% node classification forms a sidon set {0, 1 ,2, 4, 8}
% 1 - outside the fracture intersection
% 2 - on the the fracture interesection line outside facture nodes
% 4 - on the the fracture interesection line bounded by fracture nodes
% 8 - bellow fracture intersection
% 16 - above fracture intersection
%% this routine could probably simplified using the parametric equation

%% intializing node
node_classification = zeros(size(coords,1),1);
%%  predefining parameters
k = 1000;
tol = 10^-13;
teta = pi/2;
R = [cos(teta) -sin(teta); sin(teta) cos(teta)];
%% computing edges above and bellow the extended fracture intersection line 
m = (p2(:,2) - p1(:,2))/(p2(:,1) - p1(:,1));
b = p1(:,2) - m*p1(:,1);
yc = coords(:,1)*m + b;
sol = (yc - coords(:,2));
node_classification(sol > 0) = 16;
node_classification(sol < 0) = 8;

%% checking  which nodes are not inside the polygon formed by
% p1 + normal, p2 + normal, p2 - normal, p1 - normal
vec = (p2-p1)';
vec = vec/norm(vec);
normal = R*vec;
k_norm = normal'*k;
check_poly = [p1+k_norm; p2 + k_norm; p2 - k_norm; p1 - k_norm];
outside = ...
    inpolygon(coords(:,1),coords(:,2),check_poly(:,1),check_poly(:,2));
node_classification(outside==0) = 1; 

%% checking nodes bounded by fracture intersection line
p1c = repmat(p1,size(coords,1),1);
p2c = repmat(p2,size(coords,1),1);
t = (coords - p1c)./(p2c-p1c);
in_line = abs(t(:,1) - t(:,2)) < tol;
in_line_bounded = in_line;
in_line_bounded(in_line_bounded) = (t(in_line) >= 0) & (t(in_line) <= 1);
in_line_not_bounded = in_line;
in_line_not_bounded(in_line_bounded) = false;
node_classification(in_line_bounded) = 4;
node_classification(in_line_not_bounded) = 2;
end
