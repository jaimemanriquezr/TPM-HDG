function pde_stokes_2_example(k, h, mu, rho, dt, n_it, stream_function, pressure)
% Time-dependent Stokes in the stress-velocity-vorticity form
% (IncompressibleElasticityProblem) with a manufactured solution on the
% unit square (UnitSquareMesh(h)), marched n_it steps of size dt (backward
% Euler) from the initial state (initial_state: the projected initial
% velocity). Plots the solution after every step.
arguments
    k               = 2
    h               = 0.2/3
    mu              = 1e-3
    rho             = 1e-4
    dt              = 1e-2
    n_it            = 10
    stream_function = "1/2 * x^2*y^2 * cos(pi*t)"
    pressure        = "0.0"
end

get_solutions = @IncompressibleElasticityProblem.get_manufactured_solution;
mesh       = UnitSquareMesh(h);

data0 = get_solutions(stream_function, pressure, mu, rho, 0);

problem = IncompressibleElasticityProblem(viscosity = mu, density = rho, ...
    initial_velocity = data0.u, time_step = dt, time = 0, ...
    datum_dirichlet = data0.u, source_term = data0.f);
state = initial_state(problem, mesh, k);

for n = 1:n_it
    t    = n*dt;
    data = get_solutions(stream_function, pressure, mu, rho, t);
    problem.Data.datum_dirichlet = data.u;
    problem.Data.source_term     = data.f;
    problem.Data.state           = state;

    solution = solve(problem, mesh, Order= k);

    state             = get_kinematic_state(problem, solution, state);
    problem.Data.time = t;

    problem.Label = "Incompressible Elasticity | T = " + t;
    plot(solution);
end
end
