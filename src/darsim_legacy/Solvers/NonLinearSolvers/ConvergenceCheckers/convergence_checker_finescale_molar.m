% Convergence checker for FS NLSolver molar formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_finescale_molar < convergence_checker_finescale
    properties
    end
    methods 
        function PrintTitles(obj)
            % disp(['Initial Eq. 1 Residual Norm: ', num2str(obj.FirstResidualNorm(1), '%5.5e')]);
            % disp(['Initial Eq. 2 Residual Norm: ', num2str(obj.FirstResidualNorm(2), '%5.5e')]);
            % disp('');
            % disp('           ||Residual Eq1||    ||Residual Eq2||    ||delta p||    ||delta z||');
            Neq = length(obj.FirstResidualNorm);
            fprintf('Initial Eq. %1d Residual Norm: %5.5e\n', reshape([1:Neq; obj.FirstResidualNorm'], Neq, []));
            fprintf('           ');
            fprintf('||Residual Eq%d||  ', 1:Neq);
            fprintf(' ||delta p||   ');
            fprintf('||delta z%d||   ', 1:Neq-1);
            fprintf('\n');
        end
        
    end
end