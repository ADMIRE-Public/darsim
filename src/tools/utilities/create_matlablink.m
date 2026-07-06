function [matlink] = create_matlablink(text, command, varargin)
if nargin == 2
    endl = [];
else
    endl = varargin{1};
end
link_string = '<a href="matlab:%s">%s</a>';
filled_linked = sprintf(link_string,command,text);
matlink = join([filled_linked, endl]);
end
