function dsl_darcy_qk_example(N, k)
% Darcy (mixed Poisson) with Q^k on a CartesianGrid, written in the form
% language and solved on the declarative path (the Q^k twin of
% dsl_darcy_example):
%   q + grad u = 0,  div q = f  in (0,1)^2,   u = g  on the boundary,
% with the HDG flux q_hat.n = q.n + tau (u - u_hat), tau = 1.
% Prints the L2 error of u_h; expected order h^(k+1).
arguments
    N = 16
    k = 2
end
g = Datum(@(x, y) sin(pi*x).*sin(pi*y));
f = Datum(@(x, y) 2*pi^2*sin(pi*x).*sin(pi*y));

Qk = QkBasis(k); VQk = Qk.to_vector(); Mk = LegendreBasis(k);
mesh = CartesianGrid(N).add_stabilization(1);

% (q, u) are condensed element-by-element; u_hat is the global unknown
problem = FEProblem(Spaces = [FESpace([VQk, Qk], Condense = true), FESpace(Mk)]);
problem.BoundaryConditions = WeakEssentialBoundaryConditions(g, "Dirichlet");

[q, u, u_hat] = get_trial_functions(problem.Space);
[r, w, v_hat] = get_test_functions(problem.Space);
dx = Measure(preset = "dx");
dS = Measure(preset = "dS");
n = ExteriorNormal();

% local equations and the transmission (flux continuity) equation
problem.Operator = AffineOperator( ...
      dot(q, r)*dx - u*divergence(r)*dx + u_hat*dot(r, n)*dS ...
    + divergence(q)*w*dx + u*w*dS - u_hat*w*dS ...
    + f*w*dx ...
    - (dot(q, n)*v_hat*dS + u*v_hat*dS - u_hat*v_hat*dS));

solution = solve(problem, mesh, k);
e_u = compute_l2_error(solution.Unknowns.field(2), mesh, g);
fprintf("N = %d, k = %d: L2 error (u_h) = %.3e\n", N, k, e_u);
end
