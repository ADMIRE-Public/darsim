% Rachford-Rice Flash Calculator for 3 component system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Amin MB
%TU Delft
%Created: 22 July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Rachford_Rice_3component < Rachford_Rice_flash_calculator
    properties
        K_old
    end
    methods
        function obj = Rachford_Rice_3component()
        end
        function K_Wilson(obj, Status, Thermo)
            P = Status.Properties('P_2').Value;
            T = Status.T;
            N = length(P);

            comp = Thermo.comp; % A cell array containing component data
            nc = length(comp);

            % Initial guess for K-value based on Wilson correlation and
            % previous timestep
            K = zeros(N, nc);
            K_Wilson = @(Tc, Pc, w) exp(5.37 * (1 +  w) * (1 - Tc /  T) + log(P /  (Pc * 1E5))); 
            for i = 1:nc
                K(:,i) = K_Wilson(comp{i}.Tc, comp{i}.Pc, comp{i}.w);
            end

            obj.K_old = K;
        end
        function SinglePhase = Flash(obj, Status, Components, Phases, Thermo)

            % Define SinglePhase objects
            nc = length(Components);
            N = length(Status.Properties('z_1').Value);
            Z = zeros(N, nc);
            SinglePhase = zeros(N, 1);
            for i=1:nc
                Z(:, i) = Status.Properties(['z_', num2str(i)]).Value;
            end
            P = Status.Properties('P_2').Value;
            T = Status.T;

            comp = Thermo.comp; % A cell array containing component data

            Z(Z < 1E-10) = 1E-10;
            max_iter = 200;

            % Converting mass fractions > mole fractions
            % Here it is essential to solve Rachfor-Rice with mole
            % fractions as K-values are updating based on the fugacity
            % ratio and fugacity_ratio = K_mole * x_mole / y_mole
            % Z: mass fractions   z: mole fractions
            Mw = zeros(1, nc);
            for i = 1:nc
                Mw(i) = comp{i}.Mw;
            end
            z = Z ./ Mw ./ sum(Z ./ Mw, 2);

            % % Initial guess for K-value based on Wilson correlation and
            % % previous timestep
            % K = zeros(N, nc);
            % K_Wilson = @(Tc, Pc, w) exp(5.37 * (1 +  w) * (1 - Tc /  T) + log(P /  (Pc * 1E5))); 
            % for i = 1:nc
            %     K(:,i) = K_Wilson(comp{i}.Tc, comp{i}.Pc, comp{i}.w);
            % end
            % 
            % % For two-phase regions use K-value from previous timestep
            % TwoPhase = Status.Properties('SinglePhase').Value == 0;
            % K(TwoPhase, :) = obj.K_old(TwoPhase, :);
            K = obj.K_old;

            it_flash = 1; converged_flash = false; V_old = zeros(N, 1);
            while ~converged_flash
            
                F  = @(X) sum((K - 1) .* z ./ (1 + (K - 1) .* X), 2);
                dF = @(X) sum(-z .* (K - 1) .^ 2 ./(1 + X .* (K - 1)) .^ 2 , 2);
            
                V_min = zeros(N, 1);
                V_max = zeros(N, 1);
            
                index_1 = all([any(K > 1, 2), any(K < 1, 2)], 2);
                V_min(index_1) = 1 ./ (1 - max(K(index_1, :), [], 2));
                V_max(index_1) = 1 ./ (1 - min(K(index_1, :), [], 2));
                index_2 = ~index_1;
                V_max(index_2) = 1 ./ (1 - max(K(index_2, :), [], 2));
                V_min(index_2) = 1 ./ (1 - min(K(index_2, :), [], 2));
            
                F_min = F(V_min); 
                F_max = F(V_max);
                V_bs = (V_min + V_max) / 2;  % V_mid by bisection
                F_mid = F(V_bs);
                V = V_bs;
            
                % Using a hybrid algorithm of bisection and Newton-Raphson
                % Ref: Rachford–Rice equation–A look from a different perspective
                % R. Monroy-Loperena (2023)
                it_NewtonRaphson = 1;
                while ~all(abs(V_max - V_min) < 1e-15 | abs(F_mid) < 1E-30)
                    
                    % fprintf('#iter bisection-NewtonRaphson: %3d\n', it_RachfordRice);
                    V_bs = (V_min + V_max) / 2;  % Next iter solution by bisection
                    V_NR = V - F(V) ./ dF(V);    % Next iter solution by Newton-Raphson
                    
                    % Checking if the NewtonRaphson soultion is inside the window [V_min, V_max]
                    index_1 = all([V_min < V_NR, V_NR < V_max], 2); 
                    V(index_1) = V_NR(index_1);
                    index_2 = ~index_1;
                    V(index_2) = V_bs(index_2);

                    % Update F_mid
                    F_mid = F(V);
                    
                    index_N = F_min .* F_mid < 0;
                        V_max(index_N) = V(index_N);
                        F_max(index_N) = F_mid(index_N);
            
                    index_P = F_min .* F_mid > 0;
                        V_min(index_P) = V(index_P);
                        F_min(index_P) = F_mid(index_P);

                    index_0 = F_min .* F_mid == 0;
                        V_max(index_0) = V(index_0);
                        F_max(index_0) = F_mid(index_0);
                        V_min(index_0) = V(index_0);
                        F_min(index_0) = F_mid(index_0);

            
                    it_NewtonRaphson = it_NewtonRaphson + 1;
                    if it_NewtonRaphson > max_iter
                        error('Rachford-Rice/bisection-NewtonRaphson did not converged.');
                    end
                end
            
                dV = V - V_old; V_old = V;
                x = z ./ (1 + (K - 1) .* V);
                y = K .* x;
            
                fugacity_ratio = zeros(N, nc);
                for i = 1:nc
                    fugacity_ratio(:, i) = comp{i}.K_mole(P, y(:, 1)) .* x(:, i) ./ y (:, i);
                end
            
                err = sum((fugacity_ratio - 1) .^ 2, 2);
                % fprintf('Error: %10.15e\n', norm(err, inf));
            
                if all(any([err < 1e-12, all([abs(dV) < 1e-6, any([V > 1.0001, V < -0.0001], 2)], 2)], 2), 1)
                    converged_flash = true;
                    fprintf('Flash converged in %d iterations. Error: %10.15e\n', it_flash, norm(err, inf));
                    obj.K_old = K;
                else
                    K =  K .* fugacity_ratio;
                end
                it_flash = it_flash + 1;
                if it_flash > max_iter
                    keyboard;
                    % Use this to visualize F(V)
                    % i = # cell
                    % syms X; figure; fplot(sum((K(i, :) - 1) .* z(i, :) ./ (1 + (K(i, :) - 1) * X))); hold on; fplot(sym(0));
                    error('Rachford-Rice did not converged.');
                end
            
            end

            SinglePhase(V >= 1) = 1;
            V(SinglePhase == 1) = 1;
            x(SinglePhase == 1, :) = 0;
            y(SinglePhase == 1, :) = z(SinglePhase == 1, :);

            SinglePhase(V <= 0) = 2;
            V(SinglePhase == 2) = 0;
            x(SinglePhase == 2, :) = z(SinglePhase == 2, :);
            y(SinglePhase == 2, :) = 0;

            % Converting mole fractions > mass fractions
            X = x .* Mw ./ sum(x .* Mw, 2);
            Y = y .* Mw ./ sum(y .* Mw, 2);
            V = (Z - X) ./ (Y - X);
            V = mean(V, 2);

            V(SinglePhase == 1) = 1;
            X(SinglePhase == 1, :) = 0;
            Y(SinglePhase == 1, :) = Z(SinglePhase == 1, :);

            V(SinglePhase == 2) = 0;
            X(SinglePhase == 2, :) = Z(SinglePhase == 2, :);
            Y(SinglePhase == 2, :) = 0;

            % Copy it into Status object
            % Gas phase
            NI = Status.Properties('ni_1');
            NI.Value = V;
            for i = 1:nc
                X_ = Status.Properties(['x_', num2str(i),'ph1']);
                X_.Value = Y(:, i);
            end

            % Liquid phase
            NI = Status.Properties('ni_2');
            NI.Value = 1 - V;
            for i = 1:nc
                X_ = Status.Properties(['x_', num2str(i),'ph2']);
                X_.Value = X(:, i);
            end

        end
    end
end