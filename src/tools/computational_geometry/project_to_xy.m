function R = project_to_xy(points,normal)
% Takes a nx3 matrix of 3D points and returns a rotation matrix
% that projects all points onto the XY plane 
% Compute the centroid of the points
centroid = mean(points);
% Translate the points to the origin
points = points - centroid;

% Compute the normal vector to the plane of the points
%normal = null(points);

% Compute the rotation matrix that maps the normal vector to the Z-axis
v = normal(:,1);
theta = -atan2(v(1), v(2));
R1 = [cos(theta) -sin(theta) 0;      sin(theta) cos(theta)  0;      0         0           1];
v = R1*v;
phi = atan2(v(3), sqrt(v(1)^2 + v(2)^2));
R2 = [cos(phi)  0 sin(phi);      0         1 0;      -sin(phi) 0 cos(phi)];
R = R2*R1;
end