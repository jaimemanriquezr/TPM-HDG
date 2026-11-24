function r = newton_residual(pde, mesh, k, X, x_hat, indices, options)
% Residual of the hyperelastic HDG system at (X, û).
%
%   r = newton_residual(pde, mesh, k, X, x_hat)
%   r = newton_residual(pde, mesh, k, X, x_hat, indices)
%   r = newton_residual(..., nonlinear = false)    (the linear residual, m = 0)
%   r = newton_residual(..., local_rows = true)   (the local rows only,
%                       one column per element)
%
% X are the local dofs [σ; u; γ] (one column per element), x_hat the
% global hybrid vector (its Dirichlet dofs are not checked). r stacks
%   the local rows    A1 X − A2 û − Af − [m(L_h); 0; m_η(L_h)]   (every element)
%   the trace rows    Σ_K (A3 X + A4 û) − g_N                    (free dofs)
% with the linear operators of LinearElasticityProblem.get_matrices
% (inertial terms and Data.state included) and the constitutive and
% symmetry moments m of the stress remainder (see get_matrices); it
% vanishes at a solution of the nonlinear scheme. The local rows serve a
% Newton loop whose trace rows are assembled elsewhere (e.g. in an
% interface coupling); of them, the constitutive and symmetry rows are
% the only ones that hold the nonlinearity. indices is a struct with
% rhs_rows (get_assemble_indices), computed if absent.
arguments
    pde     (1,1) HyperelasticProblem
    mesh          SimplicialGrid
    k       (1,1) double
    X             double
    x_hat   (:,1) double
    indices       struct = struct.empty
    options.nonlinear (1,1) logical = true
    options.local_rows (1,1) logical = false
end
pde.Order = k;
if isempty(indices)
    indices = struct();
    [~, ~, indices.rhs_rows] = get_assemble_indices(pde, mesh, k);
end
N_ELEMENTS = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
pde.Data.iterate = [];
operators = linear_operators(pde, mesh, k);
[A1, A2, Af] = deal(operators.local_operators{:});
d2peers = size(A1, 1) - 3*d2;
ix_u = d2peers + (1:2*d2);
ix_c = [1:d2peers, d2peers + 2*d2 + (1:d2)];

rows = reshape(indices.rhs_rows, [], N_ELEMENTS);
u_K = reshape(x_hat(rows), [], 1, N_ELEMENTS);
X = reshape(X, [], 1, N_ELEMENTS);
l = pagemtimes(operators.R_u, X(ix_u, :, :)) + pagemtimes(operators.R_hat, u_K);
if options.nonlinear
    m = nonlinear_terms(pde, mesh, k, l, operators.viscous);
else
    m = 0;
end

r_local = pagemtimes(A1, X) - pagemtimes(A2, u_K) - Af;
r_local(ix_c, :, :) = r_local(ix_c, :, :) - m;
if options.local_rows
    r = reshape(r_local, [], N_ELEMENTS);
    return
end

r_K = pagemtimes(operators.A3, X) + pagemtimes(operators.A4, u_K);
[~, dir_dof, g_hat] = get_boundary_conditions(pde, mesh);
r_trace = accumarray(indices.rhs_rows(:), r_K(:), [numel(x_hat), 1]) - g_hat;
r_trace(dir_dof) = [];
r = [r_local(:); r_trace];
end
