function rotated_plane = rotate_plane_to_z(N, original_plane)
    % Compute the axis of rotation and angle of rotation for the first rotation
    axis1 = cross(N, [0 0 1]);
    cos_theta1 = dot(N, [0 0 1]);
    sin_theta1 = norm(axis1);
    axis1 = axis1 / sin_theta1;
    theta1 = atan2(sin_theta1, cos_theta1);

    % Construct the rotation matrix for the first rotation
    R1 = rotation_matrix(axis1, theta1);

    % Apply the first rotation to the original plane
    rotated_plane = R1 * original_plane;

    % Compute the axis of rotation and angle of rotation for the second rotation
    axis2 = cross(rotated_plane(:, 3), [0 0 1]);
    cos_theta2 = dot(rotated_plane(:, 3), [0 0 1]);
    sin_theta2 = norm(axis2);
    axis2 = axis2 / sin_theta2;
    theta2 = atan2(sin_theta2, cos_theta2);

    % Construct the rotation matrix for the second rotation
    R2 = rotation_matrix(axis2, theta2);

    % Apply the second rotation to the rotated plane
    rotated_plane = R2 * rotated_plane;
end

function R = rotation_matrix(axis, theta)
    % Construct the rotation matrix for a rotation around the given axis by the given angle
    c = cos(theta);
    s = sin(theta);
    C = 1 - c;
    x = axis(1);
    y = axis(2);
    z = axis(3);

    R = [c + x^2*C     x*y*C - z*s   x*z*C + y*s;
         x*y*C + z*s   c + y^2*C     y*z*C - x*s;
         x*z*C - y*s   y*z*C + x*s   c + z^2*C];
end
