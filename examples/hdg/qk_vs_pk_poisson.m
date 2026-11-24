function [results, fig] = qk_vs_pk_poisson(N, k)
% HDG for the Poisson problem -lap u = f on the unit square, u = g on its
% boundary, with exact solution u = sin(pi x) sin(pi y), solved
%   - with Q^k on the N x N CartesianGrid (hdg_poisson_cartesian), and
%   - with P^k on StructuredSimplicialGrid(N, N) (DarcyProblem),
% at the same h and tau = 1. Prints the trace unknowns and the L2 errors of
% u_h, and plots u_h for both and their difference |u_h^Q - u_h^P|.
arguments
    N = 8
    k = 2
end
u_ex = @(x, y) sin(pi*x).*sin(pi*y);
f = @(x, y) 2*pi^2*sin(pi*x).*sin(pi*y);

% Q^k
cgrid = CartesianGrid(N);
sol_q = hdg_poisson_cartesian(cgrid, k, f, u_ex);
e_q = QkBasis(k).compute_l2_error(sol_q.u, cgrid, exact_solution = u_ex);

% P^k (Darcy with identity diffusion: velocity = -grad u, pressure = u)
data = DarcyProblem.get_manufactured_solution("sin(pi*x)*sin(pi*y)", eye(2));
problem = DarcyProblem(diffusion = data.kappa, datum_dirichlet = data.p, source = data.f);
tri = StructuredSimplicialGrid(N, N);
sol_p = solve(problem, tri, Order = k);
e = sol_p.compute_l2_error({data.u, data.p}, {data.p});
e_p = e(2);
u_p = sol_p.Unknowns.field(2).dof;                       % Dubiner dof, d x n_triangles
n_trace_p = (k + 1)*(size(tri.faces, 2) - numel(tri.dirfaces));

fprintf("k = %d, N = %d | Q^k: trace dofs = %d, |u - u_h| = %.2e | P^k: trace dofs = %d, |u - u_h| = %.2e\n", ...
    k, N, sol_q.n_trace_dof, e_q, n_trace_p, e_p);
results = struct(k = k, N = N, trace_dof_Q = sol_q.n_trace_dof, error_Q = e_q, ...
    trace_dof_P = n_trace_p, error_P = e_p);

% ---- sample both solutions on a sub-triangulation of every square ----------
% Sub-triangles are cut along the anti-diagonal, so none crosses the
% diagonal of the "left" triangulation; vertices are nudged into their
% sub-triangle so that each one is located in a single triangle.
n_sub = 6;
[V, T] = sub_triangulation(n_sub);
rs = V(:, T(:));                                          % 2 x 3 n_T, unshared vertices
c = repelem(reshape(mean(reshape(rs, 2, 3, []), 2), 2, []), 1, 3);
rs = rs + 1e-9*(c - rs);
X = cgrid.get_physical_points(rs);                         % 2 x n_pts x n_el
z_q = QkBasis(k).evaluate(rs) * sol_q.u;                  % n_pts x n_el
P = reshape(X, 2, []);
[rs_tri, el] = tri.get_2d_local(P);
assert(~any(isnan(el)), "sample point not located in the triangulation")
z_p = reshape(sum(DubinerBasis(k).evaluate(rs_tri) .* u_p(:, el).', 2), size(z_q));

faces = reshape(1:numel(z_q), 3, []).';
fig = figure(Position = [100 100 1300 380]);
tl = tiledlayout(fig, 1, 3, TileSpacing = "compact");
titles = [sprintf("Q^%d on %d\\times%d squares", k, N, N), ...
          sprintf("P^%d on %d triangles", k, 2*N^2), "|u_h^Q - u_h^P|"];
values = {z_q, z_p, abs(z_q - z_p)};
for i = 1:3
    ax = nexttile(tl);
    patch(ax, Faces = faces, Vertices = [P.', values{i}(:)], ...
        FaceVertexCData = values{i}(:), FaceColor = "interp", EdgeColor = "none");
    view(ax, 3); axis(ax, "tight"); grid(ax, "on");
    colorbar(ax); title(ax, titles(i)); xlabel(ax, "x"); ylabel(ax, "y");
end
% the difference is discontinuous across element edges: show it from above,
% with the triangle edges (which contain the square edges) on top
view(ax, 2); axis(ax, "equal", "tight");
patch(ax, Faces = tri.elements.', Vertices = [tri.coordinates.', ...
    max(values{3}(:))*ones(size(tri.coordinates, 2), 1)], FaceColor = "none", ...
    EdgeColor = [1 1 1], EdgeAlpha = 0.5, LineWidth = 0.5);
title(tl, sprintf("HDG Poisson, u = sin(\\pix) sin(\\piy): |u - u_h| = %.2e (Q^%d), %.2e (P^%d)", ...
    e_q, k, e_p, k));
end

function [V, T] = sub_triangulation(n)
% (n+1)^2 lattice on [0, 1]^2, each cell cut along its anti-diagonal
[r, s] = ndgrid(linspace(0, 1, n + 1));
V = [r(:).'; s(:).'];
id = reshape(1:(n+1)^2, n + 1, n + 1);
sw = id(1:n, 1:n); se = id(2:n+1, 1:n); ne = id(2:n+1, 2:n+1); nw = id(1:n, 2:n+1);
T = [sw(:), se(:), nw(:); se(:), ne(:), nw(:)].';
end
