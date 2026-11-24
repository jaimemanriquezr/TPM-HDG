function pde_stokes_2_convergence(k, max_iterations, stream_function, pressure, viscosity)
% Stokes in the stress-velocity-vorticity form (IncompressibleElasticityProblem)
% with a manufactured solution on structured meshes of the unit square.
% Prints L2 errors of sigma_h, u_h, gamma_h, u_hat and EOC ~ k+1.
arguments
    k               = 3
    max_iterations  = 4
    stream_function = "sin(pi*x)*sin(pi*y)"
    pressure        = "x*y"
    viscosity       = 1e-3
end

data    = IncompressibleElasticityProblem.get_manufactured_solution(stream_function, pressure, viscosity);
problem = IncompressibleElasticityProblem(viscosity = data.viscosity, source_term = data.f, datum_dirichlet = data.u);

[error_sigma, error_u, error_gamma, error_uhat, h_list, N_list] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    N    = 2^i;
    mesh = StructuredSimplicialGrid(N, N);

    solution = solve(problem, mesh, Order= k);
    [error, error_hat] = solution.compute_l2_error({data.sigma, data.u, data.gamma}, {data.u});

    N_list(i)       = N;
    h_list(i)       = sqrt(2) / N;
    error_sigma(i)  = error(1);
    error_u(i)      = error(2);
    error_gamma(i)  = error(3);
    error_uhat(i)   = error_hat(1);
end

r_sigma = [NaN; compute_eoc(error_sigma, h_list)];
r_u     = [NaN; compute_eoc(error_u,     h_list)];
r_gamma = [NaN; compute_eoc(error_gamma, h_list)];
r_uhat  = [NaN; compute_eoc(error_uhat,  h_list)];

hdr = "   N  |     h    |  e(σ_h)   | r(σ_h) |   e(u_h)  | r(u_h) |  e(γ_h)   | r(γ_h) |  e(û_h)   | r(û_h)";
fprintf("\n%s\n%s\n", hdr, repmat('-', 1, strlength(hdr)));
for i = 1:max_iterations
    if i == 1
        fprintf(" %4d | %8.2e | %9.3e |   -    | %9.3e |   -    | %9.3e |   -    | %9.3e |   -\n", ...
            N_list(i), h_list(i), error_sigma(i), error_u(i), error_gamma(i), error_uhat(i));
    else
        fprintf(" %4d | %8.2e | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f\n", ...
            N_list(i), h_list(i), error_sigma(i), r_sigma(i), error_u(i), r_u(i), error_gamma(i), r_gamma(i), error_uhat(i), r_uhat(i));
    end
end
end
