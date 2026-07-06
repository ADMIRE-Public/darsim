%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outputArg1,outputArg2] = semi_plane_check_2d(p1,p2, coords)
m = (p2(:,2) - p1(:,2))/(p2(:,1) - p1(:,1));
b = p1(:,2) - m*p1(:,1);
yc = coords(:,1)*m + b;
sol = (yc - coords(:,2));

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end
