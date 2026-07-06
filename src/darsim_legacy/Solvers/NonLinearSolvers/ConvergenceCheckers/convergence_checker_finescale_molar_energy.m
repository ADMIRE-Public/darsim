% Convergence checker for FS NLSolver molar formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_finescale_molar_energy < convergence_checker_finescale
    properties
    end
    methods 
        function PrintTitles(obj)
            disp(['Initial Eq. 1 Residual Norm: ', num2str(obj.FirstResidualNorm(1), '%5.5e')]);
            disp(['Initial Eq. 2 Residual Norm: ', num2str(obj.FirstResidualNorm(2), '%5.5e')]);
            disp(['Initial Eq. 3 Residual Norm: ', num2str(obj.FirstResidualNorm(3), '%5.5e')]);
            disp('');
            disp('           ||Residual Eq1||    ||Residual Eq2||    ||Residual Eq3||    ||delta p||    ||delta z||    ||delta T||');
        end

        function converged = Check(obj, iter, residual, RHS, delta, Formulation, DiscretizationModel, State, LinearSolver)
            Nt = DiscretizationModel.N;
            
            % Compute Norms
            [ obj.ResidualNorm(iter,:), obj.RHSNorm(iter,:) ] = obj.NormCalculator.CalculateResidualNorm(residual, RHS, Nt, Formulation);
            [dp, dz, dT] = obj.NormCalculator.CalculateSolutionNorm(delta, DiscretizationModel.N, State);
            
            obj.ResidualNorm( imag(obj.ResidualNorm) ~= 0 ) = NaN;
            obj.RHSNorm     ( imag(obj.RHSNorm     ) ~= 0 ) = NaN;
            dp              ( imag(dp              ) ~= 0 ) = NaN;
            dz              ( imag(dz              ) ~= 0 ) = NaN;
            dT              ( imag(dT              ) ~= 0 ) = NaN;
            
            converged = obj.CheckConvergenceCondition(iter,dp,dz,dT);
        end
        function converged = CheckConvergenceCondition(obj,iter,dp,dz,dT)
            % disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(obj.ResidualNorm(iter,1), '%5.5e'), '      ' ...
            %                                              , num2str(obj.ResidualNorm(iter,2), '%5.5e'), '      ' ...
            %                                              , num2str(dp, '%5.5e'), '    ', num2str(dS, '%5.5e')]);
            fprintf('Iter %02d-->   ', iter);
            fprintf('%11.5e         ', obj.ResidualNorm(iter, 1:end-1));
            fprintf('%11.5e       ', obj.ResidualNorm(iter, end));
            fprintf('%11.5e    ', dp);
            fprintf('%11.5e    ', dz);
            fprintf('%11.5e   ', dT);
            fprintf('\n');


            % Initialize
            converged = 0;
            
            % check if it is stagnating
            stagnating = obj.Stagnating(obj.ResidualNorm(iter,:)'./obj.FirstResidualNorm);
            
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
                && all(dz < obj.SolutionTol(2:end-1)) ...
                && dT < obj.SolutionTol(end)

                    converged = 1;

            elseif    iter > 5 ...
                   && dp < obj.SolutionTol(1) ...
                   && all(dz < obj.SolutionTol(2:end-1)) ...
                   && all(dT < obj.SolutionTol(end)) ...
                   && all(abs(obj.ResidualNorm(end-1, :) - obj.ResidualNorm(end-2, :)) < 1e-3 * obj.ResidualNorm(end-2, :)) ...
                   && all(abs(obj.ResidualNorm(end  , :) - obj.ResidualNorm(end-1, :)) < 1e-3 * obj.ResidualNorm(end, :))
                
                    converged = 1;

            elseif    any(isnan(obj.ResidualNorm(iter, :))) ...
                   || isnan(dp) ...
                   || any(isnan(dz)) ...
                   || isnan(dT) ...
                   || stagnating

                    converged = -1;
            end
        end

        
    end
end