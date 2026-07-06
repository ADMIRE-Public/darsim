%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [WellsInfo] = build_well_info(well_array, Ngrid)
INJ_cell = cell(1);
PROD_cell = cell(1);
num_inj = 0;
num_prod = 0;
for index = 1:length(well_array)
    well = well_array{index};
    well_translated = translate_well(well, Ngrid);
    if strcmpi(well.TYPE, "INJ")
        num_inj = num_inj + 1;
        INJ_cell{num_inj} = well_translated;        
    elseif strcmpi(well.TYPE, "PROD")
        num_prod = num_prod + 1;
        PROD_cell{num_prod} = well_translated;        
    else
        error("Well Type not found.")
    end
end
WellsInfo = struct("NofWell",num_inj + num_prod, ...
    "NofInj",num_inj,"NofProd",num_prod);
%, "Inj", INJ_cell, ...
%"Prod", PROD_cell);
WellsInfo.Inj = INJ_cell';
WellsInfo.Prod = PROD_cell';
end

function well_translated = translate_well(well, Ngrid)
well_translated = struct();
NX = Ngrid(1);
NY = Ngrid(2);
NZ = Ngrid(3);
switch well.COORDINATE{1, 1}
    case('IJK_LIST')
        coord_line1 = eval(well.COORDINATE{2});
        coord_line2 = eval(well.COORDINATE{3});
        Coordinate = struct("Type", well.COORDINATE{1},"Value", [coord_line1; ...
            coord_line2]);

    case('XYZ_LIST')
        coord_line1 = eval(well.COORDINATE{2});
        coord_line2 = eval(well.COORDINATE{3});
        Coordinate = struct("Type", well.COORDINATE{1},"Value", [coord_line1; ...
            coord_line2]);

    case('CELL_INDEX_LIST')
        coord_line1 = eval(well.COORDINATE{2});
        Coordinate = struct("Type", well.COORDINATE{1},"Value", [coord_line1]);
end
Constraint = struct("Name",well.CONSTRAINT{1}, "Value",...
    well.CONSTRAINT{2});
Formula = struct("Type", well.FORMULA{1}, "Value", well.FORMULA{2});
BoundaryCondition = struct("Name",well.BOUNDARY_CONDITION{1}, "Value",...
    well.BOUNDARY_CONDITION{2});
well_translated.Coordinate = Coordinate;
well_translated.Constraint = Constraint;
well_translated.Formula = Formula;
well_translated.BoundaryCondition = BoundaryCondition;


% Check to see if the well has any schedules
try
    well_translated.Schedule = str2double(strsplit(well.SCHEDULE, ','));
catch 
    % just skip
end

% Check to see if the well has specific composition
try
    well_translated.Composition = well.COMPOSITION;
catch 
    % just skip
end

end
