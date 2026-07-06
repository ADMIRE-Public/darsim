% Compositional gas phase class considering Volume Shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Amin MB
%TU Delft
%Created: 23 July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_gas_phase_H2_CH4 < phase
    properties
        T
        % gridded interpolant object for data
        % zcf is the compressibility factor of gas phase
        zcf; % attn. zcf: compressibility factor; z: mole(mass) fraction of a component;
        dzcfdp;
        dzcfdy1;
    end
    methods
        function obj = Compositional_gas_phase_H2_CH4(mu, T) % N is number of matrix cells
            obj.T = T;
        end
        %% --- Density: (rho),   d(rho)/d(p),   d(rho)/dz ---
        function rho = ComputeDensity(obj, p, Components, Argin3)
            [x, Mw_comp, SinglePhase] = extractProperties(obj, Argin3, Components);
            [Mw_mix, y1_molar] = obj.calculateMixtureProps(x, Mw_comp);

            rho = (Mw_mix) ...
                     ./ (obj.zcf(p, y1_molar) .* 8314.462 .* obj.T ./ p );

            rho(rho == 0) = 1E-20;
        end
        function drhodp = ComputeDrhoDp(obj, p, Components, dxdp, Argin4)
            [x, Mw_comp, SinglePhase] = extractProperties(obj, Argin4, Components);
            [Mw_mix, y1_molar] = obj.calculateMixtureProps(x, Mw_comp);

            % dMw_mixture = - sum( dY_i / Mw_i ) / sum( Y_i / Mw_i ) ^ 2 
            %             = - sum( dY_i / Mw_i ) * Mw_mixture ^ 2 where Y_i is mass fraction
            dMw = - sum(dxdp(:, 1:2:6) ./ Mw_comp, 2) .* Mw_mix .^ 2;
            dMw(isnan(dMw)) = 0;

            Zcf = obj.zcf(p, y1_molar);
            dZcfdp = obj.dzcfdp(p, y1_molar);

            drhodp = (Mw_mix + p .* dMw) ./ (Zcf .* 8314.462 * obj.T) ...
                - (p .* Mw_mix .* dZcfdp) ./ (Zcf .^ 2 .* 8314.462 * obj.T);

            drhodp(SinglePhase == 2) = 0;
        end
        function drhodz = ComputeDrhoDz(obj, p, Components, dxdz, Argin4, ic)
            [x, Mw_comp, SinglePhase] = extractProperties(obj, Argin4, Components);
            [Mw_mix, y1_molar] = obj.calculateMixtureProps(x, Mw_comp);

            % dMw_mixture = - sum( dY_i / Mw_i ) / sum( Y_i / Mw_i ) ^ 2 
            %             = - sum( dY_i / Mw_i ) * Mw_mixture ^ 2 where Y_i is mass fraction
            dMwdz = - sum(dxdz(:, 1:2:6, ic) ./ Mw_comp, 2) .* Mw_mix .^ 2;
            dMwdz(isnan(dMwdz)) = 0;

            Zcf = obj.zcf(p, y1_molar);

            % d(zcf)/d(Z_j) = d(zcf)/d(Y1) * d(Y1)/d(Z_j) 
            % Because in H2-CH4-Brine thermo zcf is a function of y1
            dZcfdz = obj.dzcfdy1(p, y1_molar) .* dxdz(:, 1, ic);  

            drhodz = (p .* dMwdz) ./ (Zcf .* 8314.462 * obj.T) ...
                - (p .* Mw_mix .* dZcfdz) ./ (Zcf .^ 2 .* 8314.462 * obj.T);

            drhodz(SinglePhase == 2) = 0;

        end
        %% --- Viscosity: (mu),   d(mu)/d(p),   d(mu)/dz ---
        function mu = ComputeViscosity(obj, p, T, Argin3, Thermo)
            [x, Mw_comp, SinglePhase] = extractProperties(obj, Argin3, Thermo);
            [Mw_mix, y1_molar] = obj.calculateMixtureProps(x, Mw_comp);

            mu = Thermo.mu(p, y1_molar);
        end

        function dmudp = ComputeDmuDp(obj, p, T, Argin3, Thermo)
            [x, Mw_comp, SinglePhase] = extractProperties(obj, Argin3, Thermo);
            [Mw_mix, y1_molar] = obj.calculateMixtureProps(x, Mw_comp);

            dmudp = Thermo.dmudp(p, y1_molar);
        end

        function dmudz = ComputeDmuDz(obj, p, T, Argin3, Thermo, dxdz, c)
            [x, Mw_comp, SinglePhase] = extractProperties(obj, Argin3, Thermo);
            [Mw_mix, y1_molar] = obj.calculateMixtureProps(x, Mw_comp);

            dmudY1 = Thermo.dmudy1(p, y1_molar);
            dmudz = dmudY1 .* dxdz(:, 1, c);
        end

        %% --- Private Helper Methods ---
        function [x, Mw_comp, SinglePhase] = extractProperties(obj, Argin1, Argin2)
            % Uniform handling of input arguments
            if isa(Argin2, 'component') % Handle Thermo.comp vs Components array
                Mw_comp = [Argin2(1).Mw, Argin2(2).Mw, Argin2(3).Mw];
            else
                Mw_comp = [Argin2.comp{1}.Mw, Argin2.comp{2}.Mw, Argin2.comp{3}.Mw];
            end
            if isa(Argin1, 'status')
                N = size(Argin1.Properties('SinglePhase').Value, 1);
                nc = length(Mw_comp);
                x = zeros(N, 2*nc);
                for i = 1:nc
                    x(:, 2*i-1) = Argin1.Properties(['x_', num2str(i), 'ph1']).Value;
                end
                SinglePhase = Argin1.Properties('SinglePhase').Value;
            elseif isa(Argin1, 'double')
                x = Argin1;
                SinglePhase = zeros(size(x,1), 1);
            else
                error('Invalid argument type for property extraction.');
            end


        end

        function [Mw_mix, y1] = calculateMixtureProps(obj, x, Mw_comp)
            % Mw_mixture = 1 / sum( Y_i / Mw_i )
            Mw_mix = 1 ./ sum(x(:, 1:2:end) ./ Mw_comp, 2);
            Mw_mix(isinf(Mw_mix)) = 0;

            % y1_molar = (Y1/Mw1) / sum(Yi/Mwi)
            y1 = (x(:, 1) / Mw_comp(1)) ./ sum(x(:, 1:2:6) ./ Mw_comp(1:3), 2);
            y1(isnan(y1)) = 0; 
        end


    end    
end