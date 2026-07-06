function I = interval_intersection(II, JJ)
% Compute the intersection of two sets of intervals

% Ensure input matrices have same dimensions
assert(all(size(II) == size(JJ)), ...
    'Input matrices must have the same dimensions');

% Compute the intersection of the intervals
I = [max(II(:,1), JJ(:,1)), min(II(:,2), JJ(:,2))];
I(I(:,1) >= I(:,2), :) = 0; % Set empty intersections to [0,0]
end