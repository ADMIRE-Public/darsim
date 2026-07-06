% Compositional water phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mengjie Zhao
%TU Delft
%Created: 14 Nov 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_water_phase < phase
    properties
        Pb = 1.0e7;
        Pst = 1.01325e5;
        m
        q
        cw = 4e-10; %compressibility of water
    end
    methods
        function obj = Compositional_water_phase(mu)
            obj.mu = mu;
            % obj.sr = 0.2; % it's assigned from input file
            % obj.m = 0.2/(obj.Pb/obj.Pst - 1) * 1/obj.Pst;
            % obj.q = 1 - obj.m;
        end
        function rho = ComputeDensity(obj, p, Components, Argin3)
            N = size(p, 1);
            nc = size(Components, 2);
            if isa(Argin3, 'status')
                x = zeros(N, 2*nc);
                for i = 1:nc
                    x(:, 2*i) = Argin3.Properties(['x_', num2str(i), 'ph2']).Value;
                end
                SinglePhase = Argin3.Properties('SinglePhase').Value;
            elseif isa(Argin3, 'double')
                x = Argin3;
                SinglePhase = zeros(N, 1);
            else
                error('Confused with Argument 3.');
            end
            % rho_2 = 1/(x_1/rho_1 + x_2/rho_2 + x_3/rho_3)
            % Using actual density for water (last component) 
            % and apparent density for other components. 
            dnum = x(:, end) ./ (Components(nc).rho * exp(obj.cw * (p - obj.Pst)));
            for i = 1:nc-1
                dnum = dnum + x(:, 2*i) ./ Components(i).rho_n_phi;
            end
            rho = 1 ./ dnum;

            % To avoid singularity when there is no water 
            rho(any([SinglePhase == 1, x(:, end) == 0], 2)) = 1;
        end
        function drho = ComputeDrhoDp(obj, p, Components, dxdp, Argin4)
            % Note: d(rho_water)/d(p) = 0;
            N = size(p, 1);
            nc = size(Components, 2);
            if isa(Argin4, 'status')
                x = zeros(N, 2*nc);
                for i = 1:nc
                    x(:, 2*i) = Argin4.Properties(['x_', num2str(i), 'ph2']).Value;
                end
                SinglePhase = Argin4.Properties('SinglePhase').Value;
            elseif isa(Argin4, 'double')
                x = Argin4;
                SinglePhase = zeros(N, 1);
            else
                error('Confused with Argument 3.');
            end
            num = dxdp(:, end-1) ./ Components(nc).rho - x(:, end) .* (Components(nc).rho * obj.cw .* exp(obj.cw * (p - obj.Pst))) ./ (Components(nc).rho * exp(obj.cw * (p - obj.Pst))) .^ 2;
            dnum = x(:, end) ./ (Components(nc).rho * exp(obj.cw * (p - obj.Pst)));
            for i = 1:nc-1
                num = num + dxdp(:, 2*i) ./ Components(i).rho_n_phi;
                dnum = dnum + x(:, 2*i) ./ Components(i).rho_n_phi;
            end
            drho = - num ./ dnum .^ 2;

            % To avoid singularity when there is no water 
            drho(any([SinglePhase == 1, x(:, end) == 0], 2)) = 0;

            % rho_n_phi = Components(1).rho_n_phi;
            % rho_w = Components(2).rho;
            % drhow_dp = 0;
            % num = rho_n_phi .* rho_w;
            % den = (1-x1_ph2) .* rho_n_phi + x1_ph2 .* rho_w;
            % dnum = rho_n_phi .* drhow_dp;
            % dden = -rho_n_phi .* dx_dp(:,2) + dx_dp(:,2) .* rho_w + drhow_dp .* x1_ph2;
            % 
            % drho = (dnum .* den - dden .* num) ./ den.^2;
        end
        function drho = ComputeDrhoDz(obj, p, Components, dxdz, Argin4, ic)
            N = size(p, 1);
            nc = size(Components, 2);
            if isa(Argin4, 'status')
                x = zeros(N, 2*nc);
                for i = 1:nc
                    x(:, 2*i) = Argin4.Properties(['x_', num2str(i), 'ph2']).Value;
                end
                SinglePhase = Argin4.Properties('SinglePhase').Value;
            elseif isa(Argin4, 'double')
                x = Argin4;
            else
                error('Confused with Argument 3.');
            end
            num = dxdz(:, end-1, ic) ./ Components(nc).rho;
            dnum = x(:, end) ./ Components(nc).rho;
            for i = 1:nc-1
                num = num + dxdz(:, 2*i, ic) ./ Components(i).rho_n_phi;
                dnum = dnum + x(:, 2*i) ./ Components(i).rho_n_phi;
            end
            drho = - num ./ dnum .^ 2;

            % To avoid singularity when there is no water 
            drho(any([SinglePhase == 1, x(:, end) == 0], 2)) = 0;

            % rho_n_phi = Components(1).rho_n_phi;
            % rho_w = Components(2).rho;
            % 
            % num = rho_n_phi .* rho_w;
            % den = (1-x1_ph2) .* rho_n_phi + x1_ph2 .* rho_w;
            % 
            % dnum = 0;
            % dden = - rho_n_phi .* dx1_dz(:,2) + rho_w .* dx1_dz(:,2);
            % 
            % drho = (dnum .* den - dden .* num) ./ den.^2;
        end
    end    
end