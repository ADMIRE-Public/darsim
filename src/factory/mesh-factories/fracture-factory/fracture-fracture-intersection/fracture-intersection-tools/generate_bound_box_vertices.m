%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = generate_bound_box_vertices(p1, p2) 
    tag = logical([0 0 0; ...
                   0 0 1; ...
                   0 1 0; ...
                   0 1 1; ...
                   1 0 0; ...
                   1 0 1; ...
                   1 1 0;...
                   1 1 1]);
    is_double = isa(p1,"double");
    varargout = cell(size(tag,1),1);
    for index = 1:size(tag,1)
        ref = repmat(tag(index,:),size(p1,1),1);
        if is_double
            varargout{index} = (p1 .* ref) + (p2 .* ~ref);
        else
            varargout{index} = (p1 .* half(ref)) + (p2 .* half(~ref));
        end
    end
end
