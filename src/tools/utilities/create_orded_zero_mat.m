function [orded_zero_mat] = create_orded_zero_mat(ind_mat,ord_mat)
            %reshape_jagged_array_in_matrix(,)
        %% buggy, not to use
        [elements_id,groups_id] = find(ind_mat');
        auxmat = [groups_id ord_mat];
        [~,ref]  = sortrows(auxmat);
        elements_id = elements_id(ref);
        %f = reshape_jagged_array_in_matrix_linear(ind_mat, elements_id);
        orded_zero_mat = sparse_to_zero_ind_mat( ...
            reshape_jagged_array_in_matrix_linear(ind_mat, elements_id));
        orded_zero_mat = uint32(orded_zero_mat);
end