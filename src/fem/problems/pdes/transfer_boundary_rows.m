function [H, RHS, dir_dof] = transfer_boundary_rows(problem, mesh, local_operators, H, RHS, dir_dof, path_gradient)
%-------------------------------------------------------------------------
%Description: Transfer Path Method (Cockburn & Solano) essential boundary
%             conditions on the low-level path, for the velocity trace of
%             the fluid problems. For a Dirichlet face b of the mesh
%             boundary Gamma_h with adjacent element K, x on b and
%             x_bar = mapping(x) on the true boundary Gamma, m = x_bar - x,
%               u(x_bar) = u(x) + int_0^1 G(x + s m) m ds,
%             with G the reconstructed velocity gradient of the local
%             unknowns, extrapolated polynomially from K. Tested against
%             the trace basis on b this gives, per component,
%               M_b uhat_b + T_b (A1\A2)_K [uhat_K; rho_K] = G_b - T_b (A1\Af)_K
%             which replaces the strong Dirichlet rows of uhat on b; those
%             dofs become free. The hybrid layout is that of
%             get_assemble_indices: [uhat_x; uhat_y; rho]. Convection needs
%             nothing extra (beta lives inside A1, A2).
%             path_gradient(mesh, rs, K_list, m) -> (n_pts x n_local x n_b x 2)
%             is the class-specific reconstruction [G_phi(x) m]_a per local
%             basis function phi (see apply_transfer_conditions).
%-------------------------------------------------------------------------
arguments
    problem (1,1) FEProblem
    mesh SimplicialGrid
    local_operators (3,1) cell
    H (:,:) double
    RHS (:,1) double
    dir_dof (:,1) double
    path_gradient function_handle
end
data = problem.Data;
mapping = data.transfer_mapping;
if isa(mapping, "LevelSetCurve"); mapping = @(X) mapping.project(X); end
g = data.datum_dirichlet;

k = problem.Order;
d1 = k + 1;
N_FACES = size(mesh.faces, 2);
N_DOF = d1*N_FACES;
[A1, A2, Af] = deal(local_operators{:});
n_local = size(A1, 1);

bfaces = mesh.dirfaces(:).';
n_b = numel(bfaces);
if n_b == 0; return; end
[K_b, sn_b] = get_element_by_face(mesh, bfaces);
if iscell(K_b)
    K_b = cellfun(@(C) C(1), K_b);  sn_b = cellfun(@(C) C(1), sn_b);
end
K_b = K_b(:).';  sn_b = sn_b(:).';

% face quadrature (2k+2) and path quadrature (k+1 points)
[s_f, w_f] = get_1d_reference_quadrature(2*k + 2);
n_q = length(s_f);
[s_p, w_p] = get_1d_reference_quadrature(k);
n_p = length(s_p);
l = legendre_basis(k, s_f(:));                    % n_q x d1
W_l = w_f(:) .* l;

x_face = get_physical_points(mesh, s_f);          % 2 x n_q x N_FACES
x_face = x_face(:, :, bfaces);
x_bar = reshape(mapping(reshape(x_face, 2, [])), 2, n_q, n_b);
m = x_bar - x_face;
x_path = reshape(x_face, 2, 1, n_q, n_b) + reshape(s_p, 1, n_p) .* reshape(m, 2, 1, n_q, n_b);
x_path = reshape(x_path, 2, n_p*n_q, n_b);
m_path = reshape(repmat(reshape(m, 2, 1, n_q, n_b), 1, n_p), 2, n_p*n_q, n_b);
rs_path = mesh.get_2d_local(x_path, K_b(:));

Gm = path_gradient(mesh, rs_path, K_b, m_path);   % n_p*n_q x n_local x n_b x 2
Gm = reshape(Gm, n_p, n_q, n_local, n_b, 2);
PI = reshape(sum(w_p(:) .* Gm, 1), n_q, n_local, n_b, 2);   % path integrals
area_b = reshape(mesh.area(bfaces), 1, 1, n_b);
T = area_b .* pagemtimes(W_l, 'transpose', PI, 'none');    % d1 x n_local x n_b x 2
M_face = area_b .* pagemtimes(W_l, 'transpose', l, 'none'); % d1 x d1 x n_b
G_face = zeros(2*d1, 1, n_b);
for a = 1:2
    g_val = g{a}(reshape(x_bar(1, :, :), n_q, n_b), reshape(x_bar(2, :, :), n_q, n_b)) ...
          + zeros(n_q, n_b);
    G_face((a-1)*d1 + (1:d1), 1, :) = area_b .* reshape(W_l.' * g_val, d1, 1, n_b);
end

% local solver of the adjacent elements
S12 = pagemldivide(A1(:, :, K_b), A2(:, :, K_b));   % n_local x (6 d1 + 1) x n_b
S1f = pagemldivide(A1(:, :, K_b), Af(:, :, K_b));   % n_local x 1 x n_b
TT = [T(:, :, :, 1); T(:, :, :, 2)];                % 2 d1 x n_local x n_b
H_b = pagemtimes(TT, S12);                          % 2 d1 x (6 d1 + 1) x n_b
g_b = G_face - pagemtimes(TT, S1f);                 % 2 d1 x 1 x n_b
for b = 1:n_b   % + M_face on the face's own dofs
    for a = 1:2
        r = (a-1)*d1 + (1:d1);
        c = (a-1)*3*d1 + (sn_b(b) - 1)*d1 + (1:d1);
        H_b(r, c, b) = H_b(r, c, b) + M_face(:, :, b);
    end
end

% scatter: rows = the face's dofs, cols = the element's hybrid dofs
row_dof = [(bfaces - 1)*d1 + (1:d1).'; (bfaces - 1)*d1 + N_DOF + (1:d1).'];   % 2 d1 x n_b
face_K = mesh.boundaries(:, K_b);                                             % 3 x n_b
face_dof = reshape((face_K(:).' - 1)*d1 + (1:d1).', 3*d1, n_b);              % faces 1..3, d1 each
col_dof = [face_dof; face_dof + N_DOF; 2*N_DOF + K_b];                        % (6 d1 + 1) x n_b
rows = repmat(reshape(row_dof, 2*d1, 1, n_b), 1, 6*d1 + 1);
cols = repmat(reshape(col_dof, 1, 6*d1 + 1, n_b), 2*d1, 1);
H_T = sparse(rows(:), cols(:), H_b(:), size(H, 1), size(H, 2));

transfer_dof = row_dof(:);
H(transfer_dof, :) = H_T(transfer_dof, :);
RHS(transfer_dof) = g_b(:);
dir_dof = setdiff(dir_dof, transfer_dof);
end
