function [output] = create_multiple_range(low_lim,up_lim)
ranges = arrayfun(@(x, y) x:y, low_lim, up_lim, 'UniformOutput', false)';
output = cell2mat(ranges)';
end