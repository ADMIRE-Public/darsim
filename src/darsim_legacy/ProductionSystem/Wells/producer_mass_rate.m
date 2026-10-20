% Producer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef producer_mass_rate < producer
    properties
    end
    methods
        function obj = producer_mass_rate(PI, coord, rate)
            obj@producer(PI, coord)
            obj.qv = -rate;
            obj.p = 0;
        end
        function AdjustConstraint(obj, GravityModel, rhoT, h)

        end
        function UpdateState(obj, State, K, Mob, FluidModel, Time)
            obj.QPhases = zeros(size(obj.QPhases));
            obj.QComponents = zeros(size(obj.QComponents));
            obj.Qh = zeros(size(obj.Qh));
            
            % Check if the well is open/close
            isOpen = mod(sum(Time >= obj.Schedule), 2);
            if isOpen
                p = State.Properties(['P_',num2str(FluidModel.NofPhases)]);
                obj.p = obj.qv ./ (obj.PI .* sum(Mob(obj.Cells, :), 2) .* K(obj.Cells)) + p.Value(obj.Cells);
                f = Mob(obj.Cells, :) ./ sum(Mob(obj.Cells,:), 2);
                switch(FluidModel.name)
                    case{'SinglePhase','Immiscible'}
                        for i = 1:FluidModel.NofPhases
                            obj.QPhases(:,i) = f(:, i) .* obj.qv;
                        end
                    case{'Compositional','BlackOil'}
                        for i = 1:FluidModel.NofPhases
                            obj.QPhases(:,i) = f(:, i) .* obj.qv;
                        end
                        for j=1:FluidModel.NofComponents
                            for phase=1:FluidModel.NofPhases
                                x = State.Properties(['x_', num2str(j), 'ph', num2str(phase)]);
                                obj.QComponents(:, j) = obj.QComponents(:, j) + x.Value(obj.Cells) .* obj.QPhases(:,phase);
                            end
                        end
                    case{'Compositional_Energy'}
                        for i = 1:FluidModel.NofPhases
                            obj.QPhases(:,i) = f(:, i) .* obj.qv;
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
        function [dQdp, dQdS] = dQPhasesdPdS(obj, State, K, Mob, dMob, drho, NofPhases)
            p = State.Properties(['P_',num2str(NofPhases)]);
            dQdp = zeros(length(obj.Cells), NofPhases);
            dQdS = zeros(length(obj.Cells), NofPhases * (NofPhases - 1));
            for i = 1:NofPhases
                dQdp(:, i) = - Mob(obj.Cells,i) * obj.PI .* K(obj.Cells) + drho(obj.Cells, i) .* Mob(obj.Cells, i) * obj.PI .* K(obj.Cells).* (obj.p - p.Value(obj.Cells));
                dQdS(:, i) = dMob(obj.Cells, i) * obj.PI .* K(obj.Cells).* (obj.p - p.Value(obj.Cells));
            end
        end
        function [A, rhs] = AddToPressureSystem(obj, Mob, K, A, rhs)
            a = obj.Cells;
            for ii=1:length(a)
                A(a(ii),a(ii)) = A(a(ii),a(ii)) + obj.PI * K(a(ii)) .* Mob(a(ii));
                rhs(a(ii)) = rhs(a(ii)) + obj.PI * K(a(ii)) .* Mob(a(ii)) .* obj.p;
            end
        end
        function Q = TotalFlux(obj, Q, p, K, Mob)
            Q(obj.Cells) = Q(obj.Cells) + obj.PI .* K(obj.Cells) .* Mob(obj.Cells,1) .* (obj.p - p(obj.Cells));
        end
    end
end