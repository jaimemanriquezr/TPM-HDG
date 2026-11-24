function sol = hdg_poisson_cartesian(mesh, k, f, g, options)
%HDG_POISSON_CARTESIAN  HDG for the Poisson problem with Q^k on a CartesianGrid.
%
%   sol = hdg_poisson_cartesian(mesh, k, f, g, tau = 1)
%
%   Solves, in mixed form,
%       q + grad u = 0,  div q = f  in the grid's rectangle,  u = g on its boundary,
%   with q in [Q^k]^2, u in Q^k (QkBasis) on each element, the trace u_hat in
%   P^k (LegendreBasis) on each face, and the HDG flux
%   q_hat.n = q.n + tau (u - u_hat). The formulation is the one of
%   dsl_darcy_example: per element K, for all (r, w),
%       (q, r) - (u, div r) + <u_hat, r.n> = 0
%       (div q, w) + <tau (u - u_hat), w> = (f, w)
%   and on each interior face  sum_K <q.n + tau (u - u_hat), mu> = 0.
%   (q, u) are condensed element by element; u_hat on boundary faces is
%   the L2 projection of g. This is hdg_solve_cartesian with epsilon = 1,
%   beta = 0 and a scalar tau.
%
%   f, g are function handles of (x, y). sol has fields
%     u      Dimension-by-n_el coefficients of u (QkBasis(k))
%     q      2*Dimension-by-n_el coefficients of q, [q_x; q_y]
%     u_hat  (k+1)-by-n_faces coefficients of u_hat (LegendreBasis(k))
%     n_trace_dof   number of global trace unknowns (interior faces)
%
%See also: hdg_solve_cartesian, CartesianGrid, QkBasis, dsl_darcy_example
arguments
    mesh CartesianGrid
    k (1,1) double {mustBeNonnegative, mustBeInteger}
    f function_handle
    g function_handle
    options.tau (1,1) double = 1
end
sol = hdg_solve_cartesian(mesh, k, f, g, epsilon = 1, beta = [0; 0], tau = options.tau);
end
