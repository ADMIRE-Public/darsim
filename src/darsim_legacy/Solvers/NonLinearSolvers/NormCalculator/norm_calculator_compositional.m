% Norm calculator comp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_compositional < norm_calculator
    methods
        function [ResidualNorm, RHSNorm] = CalculateResidualNorm(obj, residual, RHS, N, Formulation)
            ResidualNorm = zeros(length(obj.FirstResidualNorm),1);
            RHSNorm = zeros(length(obj.FirstRHSNorm)     ,1);
            n_comp = Formulation.NofComponents;
            for i=1:n_comp
                ResidualNorm(i) = norm(residual((i-1)*N + 1:i*N), 2);
                RHSNorm(i) = norm(RHS     ((i-1)*N + 1:i*N), 2);
            end
            % Equilibrium = norm(residual(N*n_comp+1:end), inf);
        end
        function [dp, dz] = CalculateSolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N), inf)/max(State.Properties('P_2').Value);
            dz = vecnorm(reshape(delta(N+1:end), N, []), inf, 1);
            
%             dp = norm(delta(1:N), 2)/max(State.Properties('P_2').Value);
%             dz = norm(delta(N+1:end),2)/max(State.Properties('z_1').Value);
            
            %[a, b] = max(norm(dS));
            %disp(['max is ', num2str(a), ' in cell ', num2str(b)]);
        end
    end
end