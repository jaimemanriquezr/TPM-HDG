function pde_darcy_example(N, k, exact_solution, diffusion)
% Darcy with a manufactured solution and Dirichlet data on an N x N
% structured mesh of the unit square. Prints the L2 errors of u_h, p_h
% and p_hat.
arguments
    N = 20
    k = 1;
    exact_solution = "sin(pi*x)*sin(pi*y)";
    diffusion = eye(2);
end

data = DarcyProblem.get_manufactured_solution(exact_solution, diffusion);
problem = DarcyProblem(diffusion = data.kappa, datum_dirichlet = data.p, source = data.f);
mesh = StructuredSimplicialGrid(N, N);
solution = solve(problem, mesh, Order= k);
[e, e_hat] = solution.compute_l2_error({data.u, data.p}, {data.p});

fprintf("L2 error (u_h): %g | L2 error (p_h): %g | L2 error (p_hat): %g\n", e(1), e(2), e_hat);
end
