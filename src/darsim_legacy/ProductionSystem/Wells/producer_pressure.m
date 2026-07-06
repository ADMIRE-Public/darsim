% Producer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef producer_pressure < producer
    properties
    end
    methods
        function obj = producer_pressure(PI, coord, pressure)
            obj@producer(PI, coord)
            obj.p = pressure;
        end
        function AdjustConstraint(obj, GravityModel, rhoT, h)
            rho = rhoT(obj.Cells);
            obj.BHPDepth = max(h);
            obj.p = obj.p - rho .* GravityModel.g .*(obj.BHPDepth - h);
        end
        function UpdateState(obj, State, K, Mob, FluidModel, Time)
            obj.QPhases = zeros(size(obj.QPhases));
            obj.QComponents = zeros(size(obj.QComponents));
            
            % Check if the well is open/close
            isOpen = mod(sum(Time >= obj.Schedule), 2);
            if isOpen
                switch(FluidModel.name)
                    case{'SinglePhase','Immiscible'}
                        for i = 1:FluidModel.NofPhases
                            p = State.Properties('P_2');
                            rho = State.Properties(['rho_', num2str(i)]);
                            obj.QPhases(:,i) = rho.Value(obj.Cells) .* Mob(obj.Cells, i) .* obj.PI .* K(obj.Cells) .* (obj.p - p.Value(obj.Cells));
                        end
                    case{'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                        for i = 1:FluidModel.NofPhases
                            p = State.Properties('P_2');
                            rho = State.Properties(['rho_', num2str(i)]);
                            h = State.Properties(['h_', num2str(i)]);
                            obj.QPhases(:,i) = rho.Value(obj.Cells) .* Mob(obj.Cells, i) .* obj.PI .* K(obj.Cells) .* (obj.p - p.Value(obj.Cells));
                            obj.Qh(:,i) = h.Value(obj.Cells) .* obj.QPhases(:,i);
                        end
                    case{'Compositional','BlackOil'}
                        for i = 1:FluidModel.NofPhases
                            p = State.Properties('P_2');
                            rho = State.Properties(['rho_', num2str(i)]);
                            obj.QPhases(:,i) = rho.Value(obj.Cells) .* Mob(obj.Cells, i) .* obj.PI .* K(obj.Cells) .* (obj.p - p.Value(obj.Cells));
                        end
                        for j=1:FluidModel.NofComponents
                            for phase=1:FluidModel.NofPhases
                                x = State.Properties(['x_', num2str(j), 'ph', num2str(phase)]);
                                obj.QComponents(:, j) = obj.QComponents(:, j) + x.Value(obj.Cells) .* obj.QPhases(:,phase);
                            end
                        end
                    case{'Compositional_Energy'}
                        for i = 1:FluidModel.NofPhases
                            p = State.Properties('P_2');
                            rho = State.Properties(['rho_', num2str(i)]);
                            obj.QPhases(:,i) = rho.Value(obj.Cells) .* Mob(obj.Cells, i) .* obj.PI .* K(obj.Cells) .* (obj.p - p.Value(obj.Cells));
                        end
                        for j=1:FluidModel.NofComponents
                            for phase=1:FluidModel.NofPhases
                                x = State.Properties(['x_', num2str(j), 'ph', num2str(phase)]);
                                obj.QComponents(:, j) = obj.QComponents(:, j) + x.Value(obj.Cells) .* obj.QPhases(:,phase);
                            end
                        end
                        for ph = 1:FluidModel.NofPhases
                            h = State.Properties(['h_', num2str(ph)]);
                            obj.h(:,ph) = h.Value(obj.Cells);
                            obj.Qh(:,ph) = obj.h(:,ph) .* obj.QPhases(:,ph);
                        end
                    otherwise
                        error('This fluid model is not implemented for producers');
                end
            end
        end
        %%
        % Mass Flux Derivatives
        function dQdp = ComputeWellMassFluxDerivativeWithRespectToPressure(obj, State, K, Mob, drhodp, dMobdp, NofPhases) % need perforated cell properties
            % Q = rho * PI * K * Mob * (pWell - pCell);
            dQdp = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                p = State.Properties(['P_',num2str(NofPhases)]).Value;
                rho = State.Properties(['rho_',num2str(i)]).Value;
                dQdp(:, i) = rho(obj.Cells)      .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (-1                  ) + ...
                             drhodp(obj.Cells,i) .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells)) + ...
                             rho(obj.Cells)      .* obj.PI .* K(obj.Cells) .* dMobdp(obj.Cells,i) .* (obj.p - p(obj.Cells)) ;               
            end
        end
        function dQdS = ComputeWellMassFluxDerivativeWithRespectToSaturation(obj, State, K, dMob, NofPhases)
            dQdS = zeros(length(obj.Cells), NofPhases * (NofPhases - 1));
            for i = 1:NofPhases
                p = State.Properties(['P_',num2str(NofPhases)]).Value;
                rho = State.Properties(['rho_', num2str(i)]).Value;
                dQdS(:, i) = rho(obj.Cells) .* dMob(obj.Cells, i) .* obj.PI .* K(obj.Cells) .* (obj.p - p(obj.Cells));
            end
        end
        function dQdT = ComputeWellMassFluxDerivativeWithRespectToTemperature(obj, State, K, Mob, drhodT, dMobdT, NofPhases) % need perforated cell properties
            dQdT = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                p = State.Properties(['P_',num2str(NofPhases)]).Value;
                rho = State.Properties(['rho_',num2str(i)]).Value;
                dQdT(:, i) = obj.PI .* K(obj.Cells) .* (obj.p - p(obj.Cells)) .* ( rho(obj.Cells) .* dMobdT(obj.Cells,i) + Mob(obj.Cells,i) .* drhodT(obj.Cells,i));
            end
        end
        function dQdh = ComputeWellMassFluxDerivativeWithRespectToEnthalpy(obj, State, K, Mob, drhodh, dMobdh, NofPhases)
            % Q = rho * PI * K * Mob * (pWell - pCell);
            dQdh = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                p = State.Properties(['P_',num2str(NofPhases)]).Value;
                rho = State.Properties(['rho_',num2str(i)]).Value;
                dQdh(:, i) = drhodh(obj.Cells,i) .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells)) + ...
                             rho(obj.Cells)      .* obj.PI .* K(obj.Cells) .* dMobdh(obj.Cells,i) .* (obj.p - p(obj.Cells));
            end
        end
        %%
        % Heat Flux Derivatives
        function dQhdp = ComputeWellHeatFluxDerivativeWithRespectToPressure(obj, State, K, Mob, drhodp, dhpdp, dMobdp, NofPhases) % need perforated cell properties
            % Qh = rho * h * PI * K * Mob * (pWell - pCell);
            dQhdp = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                p = State.Properties(['P_',num2str(NofPhases)]).Value;
                rho = State.Properties(['rho_',num2str(i)]).Value;
                h = State.Properties(['h_',num2str(i)]).Value;
                dQhdp(:, i) = rho(obj.Cells)      .* h(obj.Cells)       .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (-1                  )  + ...
                              drhodp(obj.Cells,i) .* h(obj.Cells)       .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells))  + ...
                              rho(obj.Cells)      .* dhpdp(obj.Cells,i) .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells))  + ...
                              rho(obj.Cells)      .* h(obj.Cells)       .* obj.PI .* K(obj.Cells) .* dMobdp(obj.Cells,i) .* (obj.p - p(obj.Cells))        ;
            end
        end
        function dQhdT = ComputeWellHeatFluxDerivativeWithRespectToTemperature(obj, State, K, Mob, drhodT, dhdT, dMobdT, NofPhases) % need perforated cell properties
            p = State.Properties(['P_',num2str(NofPhases)]).Value;
            dQhdT = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                rho = State.Properties(['rho_',num2str(i)]).Value;
                h = State.Properties(['h_',num2str(i)]).Value;
                dQhdT(:, i) = obj.PI .* K(obj.Cells) .* (obj.p - p(obj.Cells)) .* ...
                    ( h(obj.Cells) .* rho(obj.Cells) .* dMobdT(obj.Cells,i) + ...
                    h(obj.Cells) .* Mob(obj.Cells,i) .* drhodT(obj.Cells,i) + ...
                    rho(obj.Cells) .*  Mob(obj.Cells,i) .* dhdT(obj.Cells,i));
            end
        end
        function dQhdh = ComputeWellHeatFluxDerivativeWithRespectToEnthalpy(obj, State, K, Mob, drhodh, dhpdh, dMobdh, NofPhases)
            % Qh = rho * h * PI * K * Mob * (pWell - pCell);
            dQhdh = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                p = State.Properties(['P_',num2str(i)]).Value;
                rho = State.Properties(['rho_',num2str(i)]).Value;
                h = State.Properties(['h_',num2str(i)]).Value;
                dQhdh(:, i) = drhodh(obj.Cells,i) .* h(obj.Cells)       .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells)) + ...
                              rho(obj.Cells)      .* h(obj.Cells)       .* obj.PI .* K(obj.Cells) .* dMobdh(obj.Cells,i) .* (obj.p - p(obj.Cells)) + ...
                              rho(obj.Cells)      .* dhpdh(obj.Cells,i) .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells)) ;         
            end
        end
        function q = TotalFlux(obj, q, p, K, Mob)
            q(obj.Cells) = q(obj.Cells) + obj.PI .* K(obj.Cells) .* Mob(obj.Cells,1) .* (obj.p - p(obj.Cells));
        end
    end
end