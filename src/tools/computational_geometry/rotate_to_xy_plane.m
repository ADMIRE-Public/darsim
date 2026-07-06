function R = rotate_to_xy_plane(vec)
% Takes a 3D vector vec and returns a 3x3 rotation matrix R
% that when vec is rotated by R, it is projected onto the XY plane.

% Normalize vec to unit length
vec = vec / norm(vec);

% Compute the normal vector of the plane spanned by vec and the z-axis
if vec(3) == 1
    % vec is already parallel to the XY plane, so return identity matrix
    R = eye(3);
else
    n = cross(vec, [0 0 1]); % cross product with z-axis
    n = n / norm(n);
    
    % Compute rotation angle using dot product
    cos_theta = dot(vec, [0 0 1]);
    sin_theta = norm(cross(vec, [0 0 1]));
    theta = atan2(sin_theta, cos_theta);
    
    % Compute rotation matrix using rotation formula
    I = eye(3);
    K = [0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
    R = I + sin(theta)*K + (1-cos(theta))*(K*K);
end
end