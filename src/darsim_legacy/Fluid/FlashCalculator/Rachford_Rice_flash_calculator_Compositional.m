% Rachford-Rice Flash Calculator base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini, Barnaby Fryer, Yuhang Wang, Mengjie Zhao
%TU Delft
%Created: 26 October 2016
%Last modified: 6 Jan 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Rachford_Rice_flash_calculator_Compositional < Kvalues_flash_calculator
    properties
        tol = 1e-5;
        
    end
    methods
        function obj = Rachford_Rice_flash_calculator_Compositional()

        end
        
        function SinglePhase = Flash(obj, Status, Components, Phases, Thermo)
            % Define SinglePhase objects
            nc = length(Components);
            nph = length(Phases);
            N = length(Status.Properties('z_1').Value);
            z = zeros(N, nc);
            SinglePhase = zeros(N, 1);
            for i=1:nc
                z(:, i) = Status.Properties(['z_', num2str(i)]).Value;
            end
            
            %% 1 Check if there are 2 components
            % x(z(:, 1) == 1, 1) = 1;  
            % x(z(:, 1) == 1, 2) = 0;
            % SinglePhase(z(:, 1) == 1) = 1; % All vapour
            % x(z(:, 1) == 0, 1) = 0;  
            % x(z(:, 1) == 0, 2) = 0;
            % SinglePhase(z(:, 1) == 0) = 2; % All liquid
            
            %% 2 Compute K-values
            k = obj.KvaluesCalculator.Compute(Status, Thermo);
            
            %% 3 Chek if we are in 2 phase region
            % 3.a: checking if it 's all liquid: checks if mix is below bubble
            BubblePointCheck = z .* k;
            BubblePointCheck = sum(BubblePointCheck, 2);         
            x(BubblePointCheck - 1 <= obj.tol, 2) = z (BubblePointCheck - 1 <= obj.tol, 1);
            x(BubblePointCheck - 1 <= obj.tol, 1) = 0; % This is to avoid having singular Jacobian matrix.
            x(BubblePointCheck - 1 <= obj.tol, 3) = 0;
            x(BubblePointCheck - 1 <= obj.tol, 4) = 1-z(BubblePointCheck - 1 <= obj.tol, 1);
            SinglePhase(BubblePointCheck - 1 <= obj.tol) = 2; % It s all liquid
            
            % 3.b: checking if it 's all vapor: checks if mix is above dew point
            DewPointCheck = z ./ k;
            DewPointCheck = sum(DewPointCheck, 2);
            x(DewPointCheck - 1 < obj.tol, 1) = z (DewPointCheck - 1 < obj.tol, 1);
            x(DewPointCheck - 1 < obj.tol, 2) = 0; % This is to avoid having singular Jacobian matrix.
            x(DewPointCheck - 1 < obj.tol, 3) = 1-z(DewPointCheck - 1 < obj.tol, 1);
            x(DewPointCheck - 1 < obj.tol, 4) = 0; 
            SinglePhase(DewPointCheck - 1 < obj.tol) = 1;  % It s all vapour
            
            %% 4. Actual Flash: solves for fv (vapor fraction)
            % obj.TwoPhase = ones(length(SinglePhase), 1);
            %TwoPhase = ones(length(SinglePhase), 1);
            %TwoPhase(SinglePhase > 0 ) = 0;
            TwoPhase = (SinglePhase == 0);
            fv = zeros(N,1);
            fv(SinglePhase == 1) = 1;   % all vapor
            fv(SinglePhase == 2) = 0;   % all liquid
            fv(SinglePhase == 0) = - ((k(SinglePhase == 0,2)-1)+ z(SinglePhase == 0,1).* (k(SinglePhase == 0,1)-k(SinglePhase == 0,2))) ./ ((k(SinglePhase == 0,1)-1) .*(k(SinglePhase == 0,2)-1));
           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %           %%%%%%%%% BISECTION %%%%%%%%%
            % find fv with bisection method
            % itLimit = 300;
            % 
            % fv = zeros(N,1);
            % fv(SinglePhase == 1) = 1;   % all vapor
            % fv(SinglePhase == 2) = 0;   % all liquid
            % idx = find(TwoPhase);
            % ntp = numel(idx);
            % if ntp > 0
            %     zt = z(idx,:);
            %     kt = k(idx,:);
            % 
            %     % Bisection brackets for two-phase cells
            %     fva = zeros(ntp,1);
            %     fvb = ones(ntp,1);
            % 
            %     for itCounter = 1:itLimit
            % 
            %         fvm = (fva + fvb)./2;
            % 
            %         Fa = sum( zt .* (kt - 1) ./ (1 + fva .* (kt - 1)), 2 );
            %         Fm = sum( zt .* (kt - 1) ./ (1 + fvm .* (kt - 1)), 2 );
            % 
            %         if max(abs(Fm)) < 1e-10
            %             fva = fvm; fvb = fvm;
            %             break;
            %         end
            % 
            %         sameSign = (Fa .* Fm) > 0;
            %         fva(sameSign) = fvm(sameSign);
            %         fvb(~sameSign) = fvm(~sameSign);
            % 
            %         if max(fvb - fva) < 1e-14
            %             disp(['Rachford-Rice converged in ', num2str(itCounter), ' iterations'])
            %             break;
            %         end
            %     end
            % 
            %     fv(idx) = 0.5*(fva + fvb);
            % end
            
            %% 5. Compute xs and ys
            % Solves for mole fractions in liquid phase
            x(TwoPhase == 1, 2) = z(TwoPhase == 1, 1) ./ (fv(TwoPhase == 1, 1) .* (k(TwoPhase == 1, 1) - 1) + 1);
            % Solves for mole fractions in gas phase
            x(TwoPhase == 1, 1) = k(TwoPhase == 1, 1) .* x(TwoPhase == 1, 2);
            x(TwoPhase == 1, 3:4) = 1 - x(TwoPhase == 1,1:2);
            
            Ni(:,1) = fv;
            Ni(:,2) = 1-fv;
            % Copy it into Status object
            for j=1:nph
                NI = Status.Properties(['ni_', num2str(j)]);
                NI.Value = Ni(:,j);
                for i=1:nc
                    X = Status.Properties(['x_', num2str(i),'ph',num2str(j)]);
                    X.Value = x(:,(i-1)*nph + j);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %           %%%%%%%%% NEWTON-RAPHSON %%%%%%%%%
%             % Initilaize variables
%             alpha = ones(N, 1);
%             fv = Status.Properties('ni_1').Value;
%             converged = 0;
%             itLimit = 10000;
%             hi = zeros(N, nc);
%             dhi = zeros(N, nc);
%             while ~converged && min(alpha) > 0.01
%                 itCounter = 0;
%                 while itCounter < itLimit && ~converged
%                     % Finds hi for each component
%                     for i=1:nc
%                         hi(:,i) = (z(:,i) .* k(:,i)) ./ (fv .* (k(:,i) - 1) + 1);
%                         % Finds the derivative of hi for each component
%                         dhi(:,i) = (z(:,i) .* (k(:,i) - 1).^2) ./ ((fv .* (k(:,i) - 1) + 1).^2);
%                     end
%                     h = sum(hi, 2) - 1;
%                     dh = - sum(dhi, 2);
%                     
%                     % Update fv
%                     h(TwoPhase == 0) = 0;
%                     dh(TwoPhase == 0) = 1;
%                     fvnew = alpha .* (-h ./ dh) + fv;
%                     
%                     fv = fvnew;
%                     if norm(h, inf) < 1e-10
%                         converged = 1;
%                         disp(['Rachford-Rice converged in ', num2str(itCounter + 1), ' iterations, with alpha ', num2str(min(alpha))])
%                     end
%                     itCounter = itCounter + 1;
%                 end
%                 alpha (abs(h) > 1e-10) = alpha (abs(h) > 1e-10)/2;
%                 fv (fv > 1) = 0.9;
%                 fv (fv < 0) = 0.1;
%                 fv (isnan(fv)) = Status.Properties('ni_1').Value(isnan(fv));
%             end
%             if ~converged
%                 [~, cellIndex] = max(abs(h));
%                 disp('Warning: Flash did not fully converge!');
%                 disp(['The residual norm of the equilibrium equation is ', num2str(norm(h, inf)), ' in cell ', num2str(cellIndex)]);
%                 fv (fv > 1) = 1;
%                 fv(fv < 0 ) = 0;
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%