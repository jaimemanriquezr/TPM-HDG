function [sigma, gamma, u_hat] = balanced_lift(obj, mesh, k, u, options)
% Stress, rotation and flux-balanced trace of a given volume displacement.
%
%   [sigma, gamma, u_hat] = balanced_lift(obj, mesh, k, u)
%   [sigma, gamma, u_hat] = balanced_lift(obj, mesh, k, u, boundary_trace = g_hat)
%
% For the volume displacement u held fixed, solves the static local
% problem of get_matrices without its momentum row:
%   c1 (σ, τ) - c2 (tr σ I, τ) + (γ, τ) = <û, τ n>_∂K - (u, ∇·τ)   ∀τ
%                                (σ, η) = 0                          ∀η
% on every element K (rows 1 and 3 of A1 [σ; u; γ] = A2 û + Af, as in
% lift_stress), together with the global flux balance on the skeleton
%   Σ_K <σ n - S(u - û), μ>_∂K = <g_N, μ>                           ∀μ
% for the trace û on the faces that are not constrained. Eliminating
% (σ, γ) element by element leaves a symmetric positive definite system
% in û, assembled with get_assemble_indices. This is the static
% condensation of the static problem with the u column moved to the
% right-hand side: for (u, û, σ, γ) solving the static problem (or a
% Newmark step, whose inertia enters only the momentum row) with the
% same boundary data, it returns the solved (σ, γ, û). No mass, history
% or source term enters. Linear in (u, boundary data).
%
% Boundary traces:
%   boundary_trace empty (default): as solve does, from
%     get_boundary_conditions: the Dirichlet faces carry the face L²
%     projection of the Dirichlet data (Data.boundary_conditions, or
%     Data.datum_dirichlet on the "Dirichlet" tag), the Neumann faces the
%     flux-balance row with the traction moments g_N, and a boundary face
%     with neither the flux-balance row with g_N = 0.
%   boundary_trace given (2*d1 x F Legendre face dofs, x block over y
%     block): û = boundary_trace on every boundary face (both components);
%     its entries on interior faces are not read.
% The problem's AdditionalConstraints are not applied.
%
%   u       volume dofs, 2*d2 x N (x block over y block)
%   sigma   PEERS stress dofs, d2peers x N (solution.Unknowns.field(1).dof)
%   gamma   rotation dofs, d2 x N
%   u_hat   skeleton dofs, 2*d1 x F (solution.Unknowns.hybrid(1).dof)
arguments
    obj   (1,1) LinearElasticityProblem
    mesh        SimplicialGrid
    k     (1,1) double
    u           double
    options.boundary_trace double = []
end
obj.Order = k;
[A1, A2, Af] = static_local_operators(obj, mesh, k);
N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
d2peers = size(A1, 1) - 3*d2;
ix_s = 1:d2peers;
ix_u = d2peers + (1:2*d2);
ix_g = d2peers + 2*d2 + (1:d2);
ix = [ix_s, ix_g];

% flux rows of get_matrices: Σ_K <σ n - S u, μ>_∂K + <S û, μ>_∂K
[~, nxLP, nyLP, tLP, tLL] = boundary_matrix(mesh, k);
O = zeros(3*d1, d2, N_ELEMENTS);
A3_s = [nxLP, O, nyLP, O, zeros(3*d1, d2peers - 4*d2, N_ELEMENTS);
        O, nxLP, O, nyLP, zeros(3*d1, d2peers - 4*d2, N_ELEMENTS)];
A3_u = -[tLP, O; O, tLP];
O = zeros(3*d1, 3*d1, N_ELEMENTS);
S_hat = [tLL, O; O, tLL];

% (σ, γ) = P û_K + q per element
u_e = reshape(u, 2*d2, 1, N_ELEMENTS);
X = pagemldivide(A1(ix, ix, :), cat(2, A2(ix, :, :), ...
    Af(ix, :, :) - pagemtimes(A1(ix, ix_u, :), u_e)));
P = X(:, 1:end-1, :);
q = X(:, end, :);
A3_sg = cat(2, A3_s, zeros(6*d1, numel(ix_g), N_ELEMENTS));
lhs_entries = pagemtimes(A3_sg, P) + S_hat;
rhs_entries = -pagemtimes(A3_sg, q) - pagemtimes(A3_u, u_e);

[lhs_rows, lhs_cols, rhs_rows] = get_assemble_indices(obj, mesh, k);
H = sparse(lhs_rows(:), lhs_cols(:), lhs_entries(:));
RHS = accumarray(rhs_rows, rhs_entries(:), [size(H, 1), 1]);

if isempty(options.boundary_trace)
    [x_hat, fixed, g_hat] = get_boundary_conditions(obj, mesh);
    x_hat = x_hat(1:size(H, 1));
    g_hat = g_hat(1:size(H, 1));
else
    assert(isequal(size(options.boundary_trace), [2*d1, N_FACES]), ...
        "LinearElasticityProblem:boundaryTraceSize", ...
        "boundary_trace must be 2*(k+1) x the number of faces.");
    N_DOF = d1*N_FACES;
    is_boundary = accumarray(mesh.boundaries(:), 1, [N_FACES, 1]) == 1;
    face_dof = reshape(((find(is_boundary) - 1)*d1 + (1:d1)).', [], 1);
    fixed = [face_dof; face_dof + N_DOF];
    trace = options.boundary_trace;
    x_hat = [reshape(trace(1:d1, :), [], 1); reshape(trace(d1+1:end, :), [], 1)];
    x_hat(setdiff(1:2*N_DOF, fixed)) = 0;
    g_hat = zeros(2*N_DOF, 1);
end
free = setdiff((1:size(H, 1)).', fixed);
RHS = RHS - H*x_hat + g_hat;
x_hat(free) = H(free, free) \ RHS(free);
u_hat = [reshape(x_hat(1:d1*N_FACES), d1, N_FACES);
         reshape(x_hat(d1*N_FACES + 1:end), d1, N_FACES)];

x = reshape(pagemtimes(P, element_traces(mesh, u_hat, d1)) + q, [], N_ELEMENTS);
sigma = x(1:d2peers, :);
gamma = x(d2peers + 1:end, :);
end
