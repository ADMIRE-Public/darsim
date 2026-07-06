%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FluidProperties = build_FluidProperties(blueprint)
FluidProperties = struct();
[FluidProperties.FluidModel, FluidProperties.NofPhases, ...
    FluidProperties.NofComponents] = get_fluid_model(...
    blueprint.GENERAL.FLUID_MODEL);
FluidProperties.mu = blueprint.FLUID.FLUID_VISCOSITY';
FluidProperties.Density = blueprint.FLUID.FLUID_DENSITY';
FluidProperties.FluidConductivity = ...
    blueprint.FLUID.FLUID_HEAT_CONDUCTIVITY';
FluidProperties.Gravity = blueprint.GENERAL.GRAVITY;
FluidProperties.Hysteresis = blueprint.GENERAL.HYSTERESIS;
FluidProperties.SpecificHeat = blueprint.FLUID.FLUID_SPECIFIC_HEAT';
FluidProperties.ComponentProperties = ...
    blueprint.COMPONENTS.COMPONENT_PROPERTIES';
FluidProperties.Comp = blueprint.FLUID.FLUID_COMPRESSIBILITY';
FluidProperties.Capillarity = ...
    generate_capillarity(blueprint.FLUID.CAPILLARITY);
FluidProperties.RelPerm = generate_RelPerm(blueprint.RELPERM_MODEL);
FluidProperties.ComponentList = blueprint.COMPONENTS.COMPONENT_LIST;
FluidProperties.Salinity = blueprint.COMPONENTS.SALINITY;
end

function [RelPerm] = generate_RelPerm(blueprint)
    RelPerm = struct();
    if size(blueprint.RELPERM, 1) == 1
        blueprint.RELPERM = {blueprint.RELPERM};
    end
    RelPerm.name = blueprint.RELPERM{1};
    if strcmp(RelPerm.name, 'Table')
        RelPerm.Files = {};
        for i = 2:length(blueprint.RELPERM)
            try 
                [~, name, ext] = fileparts(blueprint.RELPERM{i});
                if ~isempty(name) && ~isempty(ext)
                    RelPerm.Files{end+1} = blueprint.RELPERM{i};
                end
            catch
                break;
            end
        end
        
        if isempty(RelPerm.Files)
            RelPerm.Files = {"LookupTable.txt"};
            warning("No lookup table was provided for relative permeability. " + ...
                        "The default file 'LookupTable.txt' will be used.");
        else
            fprintf("\n\n");
            fprintf("----- Relative permeability lookup tables:\n");
            for i_lookup = 1:size(RelPerm.Files, 2)
                fprintf("%s\n", RelPerm.Files{i_lookup});
            end 
        end

        idx = find(cellfun(@(x) ischar(x) && strcmp(x,'srw_frac'), blueprint.RELPERM));
        if ~isempty(idx)
            RelPerm.srw_frac = blueprint.RELPERM{idx+1};
        end
    elseif strcmp(RelPerm.name, 'User-defined')
        try
            [~, ~, ext] = fileparts(blueprint.RELPERM{2});
            if strcmpi(ext, ".m")
                RelPerm.Files = blueprint.RELPERM{2};
            else
                error("Please provide a *.m MATLAB function.");
            end
        catch
            error("Please provide a *.m MATLAB function.");
        end
        i = 3;
        try
            RelPerm.s_irr = [blueprint.RELPERM{i} blueprint.RELPERM{i+1}];
            RelPerm.Kre = [blueprint.RELPERM{i+2} blueprint.RELPERM{i+3}];
            RelPerm.n = [blueprint.RELPERM{i+4} blueprint.RELPERM{i+5}];
        catch
            warning("");
        end
    else
        RelPerm.s_irr = [blueprint.RELPERM{2} blueprint.RELPERM{3}];
        RelPerm.Kre = [blueprint.RELPERM{4} blueprint.RELPERM{5}];
        RelPerm.n = [blueprint.RELPERM{6} blueprint.RELPERM{7}];
    end
end


function [capillarity] =  generate_capillarity(blueprint)
capillarity = struct();
if iscell(blueprint)
    if ~isempty(blueprint)
        capillarity.name = blueprint{1};
        capillarity.wetting = blueprint{2};

        if strcmp(capillarity.name, 'User-defined')
            try
                [~, ~, ext] = fileparts(blueprint{2});
                if strcmpi(ext, ".m")
                    capillarity.Files = blueprint{2};
                else
                    error("Please provide a *.m MATLAB function.");
                end
                try
                    capillarity.wetting = blueprint{3};
                end
            catch
                error("Please provide a *.m MATLAB function.");
            end
        end
    end
else
    capillarity.name = "disable";
end
end

function [model_type, num_phases, num_comp] = get_fluid_model(fluid_model)
model_type = fluid_model{1};
num_phases = fluid_model{2};
num_comp = fluid_model{3};
end

