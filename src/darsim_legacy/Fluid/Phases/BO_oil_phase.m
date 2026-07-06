% BO oil phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini, Yuhang Wang
%TU Delft
%Created: 14 July 2016
%Last modified: 6 Jan 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_oil_phase < phase
    properties
        Pb = 1.0e7;
        Pst = 1.01325e5;
        m
        q

        % gridded interpolant object for data
        % Rs and zcf are pressure dependent
        % Rs is the dissolution of gas component in liquid phase and is
        % implemented in BO_oil_phase
        % zcf is the compressibility factor of gas phase and is implemented
        % in BO_gas_phase
        Rs_p_table;
        dRs_p_table;
    end
    methods
        function obj = BO_oil_phase(mu)
            obj.mu = mu;
            % obj.sr = 0.2; % it's assigned from input file
            % obj.m = 0.2/(obj.Pb/obj.Pst - 1) * 1/obj.Pst;
            % obj.q = 1 - obj.m;
        end
        function rho = ComputeDensity(obj, p, Components, Rs)
            Bo_r = 1;
            % Bo_r = obj.m .* p + obj.q;
            rho = (Components(2).rho + Components(1).rho .* Rs)./Bo_r;
        end
        function drho = ComputeDrhoDp(obj, p, Components, Rs, dRs, dxdp)
            % Bo_r = obj.m .* p + obj.q;
            % dBo_r = obj.m;
            Bo_r = 1;
            dBo_r = 0;
            Num = Bo_r .* (Components(1).rho * dRs) - dBo_r .* (Components(2).rho + Components(1).rho .* Rs);
            Den = Bo_r.^2;
            drho = Num ./ Den;
        end
        function [Rs, dRs] = ComputeRs(obj, p)
            % based on pressure
            % Rs = 0.2 * p/obj.Pb + 0.2;
            % dRs = 0.2 / obj.Pb;
            
            % based on table lookup
            Rs = obj.Rs_p_table(p);
            dRs = obj.dRs_p_table(p);
            
            % Rs = 1.0e-8;
            % dRs = 0;
        end
        function drho = ComputeDrhoDz(obj, p, z, Components, SinglePhase, dxdz)
            % Bo_r = obj.m .* p + obj.q;
            Bo_r = ones(length(p), 1);
            drho = zeros(length(p), 1);
            dRs = zeros(length(p), 1);
            dRs(SinglePhase == 2) = Components(2).rho* Components(1).rho ./ ((1-z(SinglePhase == 2))*Components(1).rho).^2;
            drho(SinglePhase == 2) = Components(1).rho .* dRs((SinglePhase == 2)) ./ Bo_r (SinglePhase == 2);
        end
        function [Rs, dRs] = RsOfUnderSaturatedPhase(obj, z, Components, Rs, dRs, SinglePhase)
            Rs(SinglePhase == 2) = Components(2).rho .* z(SinglePhase == 2) ./ ((1 - z(SinglePhase == 2)) .* Components(1).rho);
            dRs(SinglePhase == 2) = 0;
        end
    end    
end