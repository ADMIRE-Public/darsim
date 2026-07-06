% Convergence checker for FS NLSolver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_finescale < convergence_checker
    properties
    end
    methods
        function PrintTitles(obj)
            disp(['Initial Eq. 1 Residual Norm: ', num2str(obj.FirstResidualNorm(1), '%5.5e')]);
            disp(['Initial Eq. 2 Residual Norm: ', num2str(obj.FirstResidualNorm(2), '%5.5e')]);
            disp('');
            disp('           ||Residual Eq1||    ||Residual Eq2||    ||delta p||    ||delta S||');
        end
        function converged = Check(obj, iter, residual, RHS, delta, Formulation, DiscretizationModel, State, LinearSolver)
            Nt = DiscretizationModel.N;
            
            % Compute Norms
            [ obj.ResidualNorm(iter,:), obj.RHSNorm(iter,:) ] = obj.NormCalculator.CalculateResidualNorm(residual, RHS, Nt, Formulation);
            [dp, dS] = obj.NormCalculator.CalculateSolutionNorm(delta, DiscretizationModel.N, State);
            
            obj.ResidualNorm( imag(obj.ResidualNorm) ~= 0 ) = NaN;
            obj.RHSNorm     ( imag(obj.RHSNorm     ) ~= 0 ) = NaN;
            dp              ( imag(dp              ) ~= 0 ) = NaN;
            dS              ( imag(dS              ) ~= 0 ) = NaN;
            
            converged = obj.CheckConvergenceCondition(iter,dp,dS,Formulation);
        end
        function converged = CheckConvergenceCondition(obj,iter,dp,dS,Formulation)
            % disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(obj.ResidualNorm(iter,1), '%5.5e'), '      ' ...
            %                                              , num2str(obj.ResidualNorm(iter,2), '%5.5e'), '      ' ...
            %                                              , num2str(dp, '%5.5e'), '    ', num2str(dS, '%5.5e')]);
            fprintf('Iter %02d-->   ', iter);
            fprintf('%11.5e       ', obj.ResidualNorm(iter, :));
            fprintf('%11.5e   ', dp, dS);
            fprintf('\n');


            % Initialize
            converged = 0;
            
            % check if it is stagnating
            stagnating = obj.Stagnating(obj.ResidualNorm(iter,:)./obj.FirstResidualNorm);

            bigDelta = obj.bigDelta(Formulation);
            
            % Check convergence
            % if ( (obj.ResidualNorm(iter,1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.RHSNorm(iter,1) < obj.ResidualTol(1)) ) && ...
            %    ( (obj.ResidualNorm(iter,2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.RHSNorm(iter,2) < obj.ResidualTol(2)) ) && ...
            %    ( dp < obj.SolutionTol(1) && dS < obj.SolutionTol(2) )
            % % if ( (obj.ResidualNorm(iter,1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.RHSNorm(iter,1) < obj.ResidualTol(1)) ) && ...
            % %    ( (obj.ResidualNorm(iter,2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.RHSNorm(iter,2) < obj.ResidualTol(2)) )
            %     converged = 1;
            % % elseif ( (obj.ResidualNorm(iter,1) < obj.ResidualTol(1)*0.5) || (obj.ResidualNorm(iter,1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)*0.5) || (obj.ResidualNorm(iter,1)/obj.RHSNorm(iter,1) < obj.ResidualTol(1)*0.5) ) && ...
            % %    ( (obj.ResidualNorm(iter,2) < obj.ResidualTol(2)*0.5) || (obj.ResidualNorm(iter,2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)*0.5) || (obj.ResidualNorm(iter,2)/obj.RHSNorm(iter,2) < obj.ResidualTol(2)*0.5) )
            % %     converged = 1;
            % elseif iter>5 && dp < obj.SolutionTol(1) && dS < obj.SolutionTol(2) && ...
            %        abs(obj.ResidualNorm(end-1,1) - obj.ResidualNorm(end-2,1)) < 1e-3 * obj.ResidualNorm(end-2,1) && ...
            %        abs(obj.ResidualNorm(end  ,1) - obj.ResidualNorm(end-1,1)) < 1e-3 * obj.ResidualNorm(end  ,1) && ...
            %        abs(obj.ResidualNorm(end-1,2) - obj.ResidualNorm(end-2,2)) < 1e-3 * obj.ResidualNorm(end-2,2) && ...
            %        abs(obj.ResidualNorm(end  ,2) - obj.ResidualNorm(end-1,2)) < 1e-3 * obj.ResidualNorm(end  ,2)
            %    converged = 1;
            % elseif (isnan(obj.ResidualNorm(iter,1)) || isnan(obj.ResidualNorm(iter,2)) || stagnating || isnan(dp) || isnan(dS))
            %     converged = -1;
            % end

            if    all(  any( [ obj.ResidualNorm(iter, :) < obj.ResidualTol; ...
                               obj.ResidualNorm(iter, :)./obj.FirstResidualNorm' < obj.ResidualTol; ...
                               obj.ResidualNorm(iter, :)./obj.RHSNorm < obj.ResidualTol ] )  ) ...
                && dp < obj.SolutionTol(1) ...
                && all(dS < obj.SolutionTol(2:end))

                    converged = 1;

            elseif    iter > 5 ...
                   && dp < obj.SolutionTol(1) ...
                   && all(dS < obj.SolutionTol(2:end)) ...
                   && all(abs(obj.ResidualNorm(end-1, :) - obj.ResidualNorm(end-2, :)) < 1e-3 * obj.ResidualNorm(end-2, :)) ...
                   && all(abs(obj.ResidualNorm(end  , :) - obj.ResidualNorm(end-1, :)) < 1e-3 * obj.ResidualNorm(end-1, :))
                
                    converged = 1;

            % elseif    iter > 5 ...
            %        && dp < obj.SolutionTol(1) ...
            %        && all(dS < obj.SolutionTol(2:end))
            % 
            %         converged = 1;

            elseif    any(isnan(obj.ResidualNorm(iter, :))) ...
                   || isnan(dp) ...
                   || any(isnan(dS)) ...
                   || stagnating

                    converged = -1;

            elseif bigDelta 
                    fprintf('##### A large Delta is detected.\n');
                    
                    converged = -1;

            end
        end
        function output = bigDelta(obj, Formulation)
            output = false;

            for i = 1:Formulation.NofComponents-1
                z_ratio = norm(Formulation.z(:,i) ./ Formulation.z_old(:,i), inf);
                output = output || (z_ratio > 1e+3) || (z_ratio < 1e-3);
                output = output && (z_ratio ~= Inf);
            end

        end
    end
end