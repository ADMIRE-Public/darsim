function R = rotate_to_z_axis_vec(N, pvec)
% Returns the rotation matrix that aligns the normal vector N with the positive Z-axis

% First, find the cross product between the normal vector and the positive Z-axis
R = zeros(size(N,1),9);
%pvec = [1 0 0];
%v = cross(N, repmat(pvec,size(N,1),1));
v = cross(N, pvec);
ref = vecnorm(v,2,2) == 0;
Imat1d = reshape(eye(3),1,[]);
R(ref,:) = repmat(Imat1d, sum(ref),1);
ref = ~ref;
if sum(ref) ~= 0
    % Then, find the rotation angle and axis using the Rodriguez formula
    term1 = dot(N(ref,:), pvec(ref,:) ,2);
    theta_vec = acosd(term1./vecnorm(N(ref,:),2,2));
    %theta = acosd(dot(N, [0; 0; 1]) / norm(N));
    u = v(ref,:) ./ vecnorm(v(ref,:),2,2);
    zcol = zeros(size(u,1),1);
    Imat1d = repmat(reshape(eye(3),1,[]), sum(ref),1);
    K = [zcol -u(:,3) u(:,2) u(:,3) zcol -u(:,1) -u(:,2) u(:,1) zcol];
    K2 = matprod1dvec(K,K);
    Rf = Imat1d + (sind(theta_vec) .* K) + ((1 - cosd(theta_vec)) .* K2);
    R(ref,:) = Rf;
end
end