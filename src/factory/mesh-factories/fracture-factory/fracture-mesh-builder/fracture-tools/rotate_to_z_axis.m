%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = rotate_to_z_axis(N)
% Returns the rotation matrix that aligns the normal vector N with the positive Z-axis

% First, find the cross product between the normal vector and the positive Z-axis
v = cross(N, [0; 0; 1]);

% If the cross product is zero, the vectors are already aligned
if norm(v) == 0
    R = eye(3);
    return
end

% Then, find the rotation angle and axis using the Rodriguez formula
theta = acosd(dot(N, [0; 0; 1]) / norm(N));
u = v / norm(v);
K = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
R = eye(3) + sind(theta) * K + (1 - cosd(theta)) * K^2;
end
