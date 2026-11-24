function r = newton_residual(pde, mesh, k, X, x_hat, indices)
% Residual of the gradient-based hyperelastic HDG system at (X, û):
% X = [H; u] local dofs (one column per element), x_hat the hybrid vector.
% r stacks the local rows (every element) and the free trace rows, with
% mesh.tau as the stabilization; it vanishes at a solution.
arguments
    pde     (1,1) HyperelasticGradientProblem
    mesh          SimplicialGrid
    k       (1,1) double
    X             double
    x_hat   (:,1) double
    indices       struct = struct.empty
end
pde.Order = k;
if isempty(indices)
    indices = struct();
    [~, ~, indices.rhs_rows] = get_assemble_indices(pde, mesh, k);
end
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
% the operators linearized at X itself: A1 X − A2 û − Af is the residual
pde.Data.iterate = struct("gradient", X(1:4*d2, :));
[local_operators, ~, ~] = get_matrices(pde, mesh, k);
[A1, A2, Af] = deal(local_operators{:});
u_K = reshape(x_hat(reshape(indices.rhs_rows, [], N_ELEMENTS)), [], 1, N_ELEMENTS);
X = reshape(X, [], 1, N_ELEMENTS);
r_local = pagemtimes(A1, X) - pagemtimes(A2, u_K) - Af;

[~, ~, ~, tLP, tLL] = boundary_matrix(mesh, k);
[~, C] = stress_terms(pde, mesh, k, reshape(X(1:4*d2, :, :), 4*d2, N_ELEMENTS));
u = X(4*d2 + 1:end, :, :);
r_K = C - [pagemtimes(tLP, u(1:d2, :, :)); pagemtimes(tLP, u(d2 + 1:end, :, :))] ...
    + [pagemtimes(tLL, u_K(1:3*d1, :, :)); pagemtimes(tLL, u_K(3*d1 + 1:end, :, :))];
[~, dir_dof, g_hat] = get_boundary_conditions(pde, mesh);
r_trace = accumarray(indices.rhs_rows(:), r_K(:), [numel(x_hat), 1]) - g_hat;
r_trace(dir_dof) = [];
r = [r_local(:); r_trace];
end
