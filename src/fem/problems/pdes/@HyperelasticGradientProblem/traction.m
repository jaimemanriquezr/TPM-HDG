function t = traction(obj, mesh, solution)
% Piola traction P̂n = P(I + H_h) n − τ (u_h − û_h) on the element faces.
%
%   t = traction(obj, mesh, solution)
%
% returns, per element (columns), the Legendre coefficients of P̂n on
% its three local faces, 6*(k+1) x N: x component (faces 1, 2, 3, k+1
% coefficients each, in the face orientation of the trace dofs) over the
% y component. mesh must carry the τ of the solve (solution.Mesh does).
arguments
    obj      (1,1) HyperelasticGradientProblem
    mesh           SimplicialGrid
    solution (1,1) FESolution
end
k = solution.Order;
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
N_ELEMENTS = size(mesh.elements, 2);
H = solution.Unknowns.field(1).dof;
u = solution.Unknowns.field(2).dof;
if iscell(u), u = cat(1, u{:}); end
u_hat = solution.Unknowns.hybrid(1).dof;
if iscell(u_hat), u_hat = cat(1, u_hat{:}); end
[~, C] = stress_terms(obj, mesh, k, H);
[~, ~, ~, tLP, tLL] = boundary_matrix(mesh, k);
u = reshape(u, 2*d2, 1, N_ELEMENTS);
u_K = element_traces(mesh, u_hat, d1);
moments = C - [pagemtimes(tLP, u(1:d2, :, :)); pagemtimes(tLP, u(d2 + 1:end, :, :))] ...
    + [pagemtimes(tLL, u_K(1:3*d1, :, :)); pagemtimes(tLL, u_K(3*d1 + 1:end, :, :))];
% the face mass matrices <L_j, L_i>_e
[~, ~, ~, ~, LL] = boundary_matrix(add_stabilization(mesh, 1), k);
O = zeros(3*d1, 3*d1, N_ELEMENTS);
t = reshape(pagemldivide([LL, O; O, LL], moments), 6*d1, N_ELEMENTS);
end
