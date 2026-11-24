function pde_elasticity_example(k, h, E, nu, K, rho, dt, n_it, u_exact)
% Time-dependent linear elasticity with a manufactured displacement on
% the unit square (UnitSquareMesh(h)), marched n_it steps of size dt from
% the initial state (initial_state: projections of the initial data).
% Plots the solution after every step.
arguments
    k       = 2
    h       = 0.2/3
    E       = 1.0
    nu      = 0.3
    K       = 0.0
    rho     = 1.1
    dt      = 1e-2
    n_it    = 10
    u_exact = ["sin(pi*x)*sin(pi*y)*cos(2*pi*t)"; "sin(pi*x)*sin(pi*y)*cos(2*pi*t)"]
end

get_solutions = @LinearElasticityProblem.get_manufactured_solution;
mesh       = UnitSquareMesh(h);

data0 = get_solutions(u_exact, E, nu, K, 0, 0);

problem = LinearElasticityProblem(young_modulus = E, poisson_ratio = nu, hooke_modulus = K, ...
    density = rho, initial_displacement = data0.u, initial_velocity = data0.u_t, time_step = dt, ...
    time = 0, datum_dirichlet = data0.u, source_term = data0.f);
state = initial_state(problem, mesh, k);

problem.Data.state = state;
solution = solve(problem, mesh, Order= k);
plot(solution);

for n = 1:n_it
    t    = n*dt;
    data = get_solutions(u_exact, E, nu, K, rho, t);
    problem.Data.datum_dirichlet = data.u;
    problem.Data.source_term     = data.f;
    problem.Data.state           = state;

    solution = solve(problem, mesh, Order= k);

    state             = get_kinematic_state(problem, solution, state);
    problem.Data.time = t;

    problem.Label = "Linear Elasticity | T = " + t;
    plot(solution);
end
end
