%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pedfm] = compute_physical_prop_pedfm(rockmesh, fracmesh, ...
    geometry,fracmeshint, edfm_geo, pedfm_geo, avg_d_rock, avg_d_frac)
    % Compute pEDFM and EDFM physical properties
    pedfm = struct();
    num_frac = size(fracmesh,1);
    %% compute infrac CI
    infrac = cell(num_frac,1);
    for index = 1:num_frac
        physical = struct();
        physical.frac = index;
        physical.K = geometry(index).Permeability;
        d1 = fracmesh(index).elements.centroid(fracmesh(index).faces. ...
            inner.upstream,:);
        d2 = fracmesh(index).elements.centroid(fracmesh(index).faces. ...
            inner.downstream,:);
        physical.d = vecnorm(d2-d1,2,2);
        p1 = fracmesh(index).vertex.coords(fracmesh(index).faces. ...
            inner.vertices(:,1),:);
        p2 = fracmesh(index).vertex.coords(fracmesh(index).faces. ...
            inner.vertices(:,2),:);
        physical.face_length = vecnorm(p2-p1,2,2);
        physical.CI = physical.face_length./physical.d;
        infrac{index} = physical;
    end
    pedfm.infrac = infrac;
    
    %% computing fracture - fracture EDFM CI transmissibilities
    fracfrac_trans = cell(num_frac,num_frac);
    for ii = 1:size(fracmeshint,1)
        for jj = ii+1:size(fracmeshint,2)
            if ~isempty(fracmeshint{ii,jj})
                physical = struct();
                local_frac = fracmeshint{ii,jj};
                local_ind1 = local_frac(:,1);
                local_ind2 = local_frac(:,2);
                int_dist = local_frac(:,3);
                auxcell = avg_d_frac{ii,jj};
                d1 = auxcell{1};
                d2 = auxcell{2};
                K1 = geometry(ii).Permeability;
                K2 = geometry(jj).Permeability;
                physical.K = (2*K1*K2)./(K1+K2);
                local_d1 = apply_filter(local_ind1, d1, unique(local_ind1));
                local_d2 = apply_filter(local_ind2, d2, unique(local_ind2));
                physical.ind1 = local_ind1;
                physical.ind2 = local_ind2;
                CI1 = 2*int_dist./local_d1;
                CI2 = 2*int_dist./local_d2;
                CIHarm = 2* (CI1.*CI2)./(CI1 + CI2);
                physical.CI =CIHarm;
                fracfrac_trans{ii,jj} = physical;
            end
        end
    end
    pedfm.fracfrac = fracfrac_trans;
    
    %% computing fracture-rock EDFM CI transmissiblities
    fracrock_edfmtrans = cell(num_frac,1);
    for ii = 1:num_frac
        physical = struct();
        vol_ind = edfm_geo.vol_ind{ii};
        frac_ind = edfm_geo.frac_ind{ii};
        shared_frac_area = edfm_geo.area{ii};
        local_d_ind = avg_d_rock{ii,1};
        local_d_value = avg_d_rock{ii,2};
        d_vol_ind = apply_filter(vol_ind,local_d_value,local_d_ind);
        CI = shared_frac_area./d_vol_ind;
        K1 = rockmesh.props.perm(vol_ind);
        K2 = geometry(ii).Permeability;
        K2 = repmat(K2,size(K1,1),1);
        L = rockmesh.element.volume(vol_ind).^(1/3);
        w = geometry(ii).Aperture;
        physical.vol_ind = vol_ind;
        physical.frac_ind = frac_ind;
        physical.K = (L + w).*(K1.*K2)./(w.*K1+L.*K2);
        physical.CI = CI;
        fracrock_edfmtrans{ii} = physical;
    end
    pedfm.fracrock_edfm = fracrock_edfmtrans;
    pedfm.pedfm_on = false;
    pedfm.alpha = zeros(size(rockmesh.face.inner.area));
    if ~isempty(fieldnames(pedfm_geo))
        %% computing fracrock pEDFM
        pedfm.pedfm_on = true;
        fracrock_pedfmtrans = cell(size(fracmeshint,1),1);
        for ii = 1:num_frac
            physical = struct();
            vol_ind = pedfm_geo.vol_ind{ii};
            frac_ind = pedfm_geo.frac_ind{ii};
            shared_frac_area = pedfm_geo.area{ii};
            local_d_ind = avg_d_rock{ii,1};
            local_d_value = avg_d_rock{ii,2};
            d_vol_ind = apply_pdefm_filter(vol_ind,local_d_ind,local_d_value);
            CI = shared_frac_area./d_vol_ind;
            K1 = rockmesh.props.perm(vol_ind);
            K2 = geometry(ii).Permeability;
            K2 = repmat(K2,size(K1,1),1);
            L = rockmesh.element.volume(vol_ind).^(1/3);
            w = geometry(ii).Aperture;
            physical.K = (L + w).*(K1.*K2)./(w.*K1+L.*K2);
            physical.CI = CI;
            physical.vol_ind = vol_ind;
            physical.frac_ind = frac_ind;
            fracrock_pedfmtrans {ii} = physical;
        end
        pedfm.fracrock_pedfm = fracrock_pedfmtrans;
    end
    %% computing alphas pEDFM
    if ~isStructEmpty(pedfm_geo)
    
        blockTol = 0.1;
        A_total = [rockmesh.face.inner.area ; rockmesh.face.bound.area];
        nFace   = numel(A_total);
    
        faceAreaLoss = accumarray( ...
            pedfm_geo.alpha_ind(:) , ...   % key
            pedfm_geo.alpha_area(:), ...   % val
            [nFace 1] , @sum , 0);         % size, func, default 0
    
        A_eff = A_total - faceAreaLoss;
        A_eff(A_eff < 0) = 0;
    
        blockedFrac = faceAreaLoss ./ A_total;
        A_eff(blockedFrac > blockTol) = 0;
    
        alpha_geom        = faceAreaLoss ./ A_total;
        alpha_geom(alpha_geom>1) = 1;
        pedfm.alpha       = alpha_geom;
        pedfm.alpha(A_total==0)=0;
    
    end

end

function [filter] = apply_filter(ind,value, all_ind)
    all_ind = double(all_ind);
    transf_matrix = sparse(all_ind, ones(length(all_ind),1), value);
    filter = full(transf_matrix(ind));   
end

function [filter] = apply_pdefm_filter(input, ind,value)
    ind = double(ind);    
    oind = ones(size(ind));
    max_value = max(ind);
    transf_matrix = sparse(ind,oind,value,max_value,1);
    filter = full(transf_matrix(input));
end
