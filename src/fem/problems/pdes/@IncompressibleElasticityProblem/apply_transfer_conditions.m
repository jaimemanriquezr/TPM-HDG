function [H, RHS, dir_dof] = apply_transfer_conditions(problem, mesh, local_operators, H, RHS, dir_dof)
%-------------------------------------------------------------------------
%Description: TPM transfer of the Dirichlet datum from the true boundary
%             to Gamma_h (transfer_boundary_rows) with the velocity
%             gradient reconstructed from the local unknowns
%             [sigma (PEERS); u; gamma; xi] as
%               grad u = c1 sigma - c2 tr(sigma) I + gamma,
%             c1 = 1/(2mu), c2 = 1/(4mu) (as trace_extension_coefficients).
%             No-op without data.transfer_mapping.
%-------------------------------------------------------------------------
if ~isfield(problem.Data, "transfer_mapping") || isempty(problem.Data.transfer_mapping)
    return
end
k = problem.Order;
mu = problem.Data.viscosity;
c1 = 1/(2*mu);  c2 = 1/(4*mu);
d2 = nchoosek(k + 2, 2);
n_sigma = 4*d2 + (k + 1)*(k > 0);
n_local = size(local_operators{1}, 1);
G_inv = get_affine_map(mesh, inverse = true);
    function Gm = path_gradient(~, rs, K_list, m)
        n_b = numel(K_list);
        n_pts = size(rs, 2);
        Gm = zeros(n_pts, n_local, n_b, 2);
        for b = 1:n_b
            r = rs(1, :, b).';  s = rs(2, :, b).';
            P = PEERS_basis(k, r, s, G_inv(:, :, K_list(b)));    % n_pts x n_sigma x 2 x 2
            A = asym_tensor_dubiner_basis(k, r, s);              % n_pts x d2 x 2 x 2
            G = c1*P;
            tr = P(:, :, 1, 1) + P(:, :, 2, 2);
            G(:, :, 1, 1) = G(:, :, 1, 1) - c2*tr;
            G(:, :, 2, 2) = G(:, :, 2, 2) - c2*tr;
            m1 = m(1, :, b).';  m2 = m(2, :, b).';
            for a = 1:2
                Gm(:, 1:n_sigma, b, a) = G(:, :, a, 1) .* m1 + G(:, :, a, 2) .* m2;
                Gm(:, n_sigma + 2*d2 + (1:d2), b, a) = A(:, :, a, 1) .* m1 + A(:, :, a, 2) .* m2;
            end
        end
    end
[H, RHS, dir_dof] = transfer_boundary_rows(problem, mesh, local_operators, H, RHS, dir_dof, @path_gradient);
end
