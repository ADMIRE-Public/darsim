function [filled_linked] = create_weblink(text, link, varargin)
if nargin == 2
    endl = newline;
else
    endl = varargin{1};
end
link_string = '<a href = "%s">%s</a>';
filled_linked = sprintf(link_string,link,text);
filled_linked = join([filled_linked, endl]);
end