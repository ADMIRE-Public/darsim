% Norm calculator comp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_compositional_energy < norm_calculator
    methods
        function [ResidualNorm, RHSNorm] = CalculateResidualNorm(obj, residual, RHS, N, Formulation)
            ResidualNorm = zeros(length(obj.FirstResidualNorm),1);
            RHSNorm = zeros(length(obj.FirstRHSNorm)     ,1);
            n_comp = Formulation.NofComponents;
            for i=1:(n_comp + 1)
                ResidualNorm(i) = norm(residual((i-1)*N + 1:i*N), 2);
                RHSNorm(i) = norm(RHS     ((i-1)*N + 1:i*N), 2);
            end
            % Equilibrium = norm(residual(N*n_comp+1:end), inf);
        end
        function [dp, dz, dT] = CalculateSolutionNorm(obj, delta, N, State)
            
            dp = norm(delta(1:N), inf)/max(State.Properties('P_2').Value);
            dz = vecnorm(reshape(delta(N+1:end-N), N, []), inf, 1);
            dT = norm(delta(end-N+1:end), inf)/max(State.Properties('T').Value);
            
            
        end
    end
end