function pde_stokes_example(N, k, stream_function, pressure, viscosity)
% Stokes in the gradient-velocity-pressure form (StokesProblem) with a
% manufactured solution on an N x N structured mesh of the unit square.
% Prints the L2 errors of L_h, u_h, p_h and u_hat.
arguments
    N               = 20
    k               = 2
    stream_function = "x^2*y^2"
    pressure        = "x*y"
    viscosity       = 1.0
end

data    = StokesProblem.get_manufactured_solution(stream_function, pressure, viscosity);
problem = StokesProblem(viscosity = data.nu, source_term = data.f, datum_dirichlet = data.u);
mesh    = StructuredSimplicialGrid(N, N);
solution = solve(problem, mesh, Order= k);
[error, error_hat] = solution.compute_l2_error({data.L, data.u, data.p}, {data.u});

fprintf("L2 error (L_h): %g | L2 error (u_h): %g | L2 error (p_h): %g | L2 error (u_hat): %g\n", ...
    error(1), error(2), error(3), error_hat(1));
end
