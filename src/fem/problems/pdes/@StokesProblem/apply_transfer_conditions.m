function [H, RHS, dir_dof] = apply_transfer_conditions(problem, mesh, local_operators, H, RHS, dir_dof)
%-------------------------------------------------------------------------
%Description: TPM transfer of the Dirichlet datum from the true boundary
%             to Gamma_h (transfer_boundary_rows) with the velocity
%             gradient L of the local unknowns [L11; L21; L12; L22; u; p; xi]
%             as the path gradient. No-op without data.transfer_mapping.
%-------------------------------------------------------------------------
if ~isfield(problem.Data, "transfer_mapping") || isempty(problem.Data.transfer_mapping)
    return
end
k = problem.Order;
d2 = nchoosek(k + 2, 2);
n_local = size(local_operators{1}, 1);
    function Gm = path_gradient(~, rs, K_list, m)
        % [L_phi m]_a: basis function phi of block L_ab contributes P_phi m_b
        n_b = numel(K_list);
        n_pts = size(rs, 2);
        Gm = zeros(n_pts, n_local, n_b, 2);
        for b = 1:n_b
            P = dubiner_basis(k, rs(1, :, b).', rs(2, :, b).');   % n_pts x d2
            m1 = m(1, :, b).';  m2 = m(2, :, b).';
            Gm(:, 0*d2 + (1:d2), b, 1) = P .* m1;    % L11 m1
            Gm(:, 1*d2 + (1:d2), b, 2) = P .* m1;    % L21 m1
            Gm(:, 2*d2 + (1:d2), b, 1) = P .* m2;    % L12 m2
            Gm(:, 3*d2 + (1:d2), b, 2) = P .* m2;    % L22 m2
        end
    end
[H, RHS, dir_dof] = transfer_boundary_rows(problem, mesh, local_operators, H, RHS, dir_dof, @path_gradient);
end
