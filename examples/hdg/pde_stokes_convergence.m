function pde_stokes_convergence(k, max_iterations, stream_function, pressure, viscosity)
% Stokes in the gradient-velocity-pressure form (StokesProblem) with a
% manufactured solution on structured meshes of the unit square. Prints
% L2 errors of L_h, u_h, p_h, u_hat and EOC ~ k+1.
arguments
    k               = 1
    max_iterations  = 8
    stream_function = "x^2*y^2"
    pressure        = "x*y"
    viscosity       = 1e-2
end

data    = StokesProblem.get_manufactured_solution(stream_function, pressure, viscosity);
problem = StokesProblem(viscosity = data.nu, source_term = data.f, datum_dirichlet = data.u);

[error_L, error_u, error_p, error_uhat, h_list, N_list] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    N    = 2^i;
    mesh = StructuredSimplicialGrid(N, N);

    solution = solve(problem, mesh, Order= k);
    [error, error_hat] = solution.compute_l2_error({data.L, data.u, data.p}, {data.u});

    N_list(i)     = N;
    h_list(i)     = sqrt(2) / N;
    error_L(i)    = error(1);
    error_u(i)    = error(2);
    error_p(i)    = error(3);
    error_uhat(i) = error_hat(1);
end

r_L    = [NaN; compute_eoc(error_L,    h_list)];
r_u    = [NaN; compute_eoc(error_u,    h_list)];
r_p    = [NaN; compute_eoc(error_p,    h_list)];
r_uhat = [NaN; compute_eoc(error_uhat, h_list)];

hdr = "   N  |     h    |   e(L_h)  | r(L_h) |   e(u_h)  | r(u_h) |   e(p_h)  | r(p_h) |  e(û_h)   | r(û_h)";
fprintf("\n%s\n%s\n", hdr, repmat('-', 1, strlength(hdr)));
for i = 1:max_iterations
    if i == 1
        fprintf(" %4d | %8.2e | %9.3e |   -    | %9.3e |   -    | %9.3e |   -    | %9.3e |   -\n", ...
            N_list(i), h_list(i), error_L(i), error_u(i), error_p(i), error_uhat(i));
    else
        fprintf(" %4d | %8.2e | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f\n", ...
            N_list(i), h_list(i), error_L(i), r_L(i), error_u(i), r_u(i), error_p(i), r_p(i), error_uhat(i), r_uhat(i));
    end
end
end
