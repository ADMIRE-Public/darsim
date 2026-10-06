% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef injector_rate < injector
    properties
     
    end
    methods
        function obj = injector_rate(PI, coord, rate, p_init, temperature, n_phases)
            obj@injector(PI, coord, n_phases)
            obj.qv = rate;
            obj.T = temperature;
            obj.p = p_init;
        end
        function AdjustConstraint(obj, GravityModel, h)
           
        end
        function UpdateState(obj, State, K, FluidModel, Time)
            obj.QPhases = zeros(size(obj.QPhases));
            obj.QComponents = zeros(size(obj.QComponents));
            obj.Qh = zeros(size(obj.Qh));
            
            % Check if the well is open/close
            isOpen = mod(sum(Time >= obj.Schedule), 2);
            if isOpen
                p = State.Properties(['P_',num2str(FluidModel.NofPhases)]); 
                %obj.p = p.Value(obj.Cells); % approximation
                obj.p = obj.qv ./ (obj.PI .* sum(obj.Mob, 2) .* K(obj.Cells)) + p.Value(obj.Cells);
                f = obj.Mob ./ sum(obj.Mob, 2);         
                switch(FluidModel.name)
                    case{'SinglePhase','Immiscible'}
                        for i = 1:FluidModel.NofPhases
                            obj.QPhases(:,i) = obj.rho(:,i) .* f(:, i) .* obj.qv;
                        end
                    case{'Compositional','BlackOil'}
                        for i = 1:FluidModel.NofPhases
                            obj.QPhases(:,i) = obj.rho(:,i) .* f(:, i) .* obj.qv;
                        end
                        for j=1:FluidModel.NofComponents
                            for phase=1:FluidModel.NofPhases
                                obj.QComponents(:, j) = obj.QComponents(:, j) + obj.x(:,(j-1)*2 + phase) .* obj.QPhases(:, phase);
                            end
                        end
                    case{'Compositional_Energy'}
                        c = obj.Composition;
                        obj.rho(:,c) = FluidModel.Thermo.rhoV_pTz_table(obj.p, obj.T, obj.z(:,c));
                        obj.Mob(:,c) = 1 ./ FluidModel.Thermo.muV_pT_table(obj.p, obj.T);
                        f = obj.Mob ./ sum(obj.Mob, 2);
                        for i = 1:FluidModel.NofPhases
                            obj.QPhases(:,i) = obj.rho(:,i) .* f(:, i) .* obj.qv;
                        end
                        for j=1:FluidModel.NofComponents
                            for phase=1:FluidModel.NofPhases
                                obj.QComponents(:, j) = obj.QComponents(:, j) + obj.x(:,(j-1)*2 + phase) .* obj.QPhases(:, phase);
                            end
                        end
                        for phase=1:FluidModel.NofPhases
                            if phase == 1
                                obj.h(:,phase) = FluidModel.Thermo.hV_pTz_table(obj.p, obj.T, obj.z(:,1));
                            else
                                obj.h(:,phase) = FluidModel.Thermo.hb_pT_table(obj.p, obj.T);
                            end
                            obj.Qh(:,phase) = obj.h(:,phase) .* obj.QPhases(:,phase);
                        end
                end
            end
        end
        function [dQdp, dQdS] = dQPhasesdPdS(obj, K, NofPhases)
            dQdp = zeros(length(obj.Cells), NofPhases);
            dQdS = zeros(length(obj.Cells), NofPhases * (NofPhases - 1));
        end
        function [A, rhs] = AddToPressureSystem(obj, K, A, rhs)
            a = obj.Cells;
            for ii=1:length(a)
                rhs(a(ii)) = rhs(a(ii)) + sum(obj.QPhases(ii,:));
            end
        end
        function q = TotalFlux(obj, q, p, K)
            q(obj.Cells) = q(obj.Cells) + obj.PI .* K(obj.Cells) .* obj.Mob(:,1) .* (obj.p - p(obj.Cells));
        end
    end
end