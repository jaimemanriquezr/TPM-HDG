function pde_darcy_convergence(k, max_iterations, exact_solution, diffusion)
% Darcy convergence study with a manufactured solution on structured
% meshes of the unit square: Neumann data on x = 0, 1, Dirichlet on
% y = 0, 1. Prints L2 errors of u_h, p_h, p_hat and EOC ~ k+1.
arguments
    k               = 1
    max_iterations  = 7;
    exact_solution  = "sin(x)*sin(y)";
    diffusion       = eye(2);
end

n_x   = @(x,y) sign(x - 1/2) .* (abs(x - 1/2) == 1/2);
n_y   = @(x,y) sign(y - 1/2) .* (abs(y - 1/2) == 1/2);
dot_n = @(u) @(x,y) u{1}(x,y).*n_x(x,y) + u{2}(x,y).*n_y(x,y);

data    = DarcyProblem.get_manufactured_solution(exact_solution, diffusion);
problem = DarcyProblem(diffusion = data.kappa, datum_dirichlet = data.p, ...
    source = data.f, datum_neumann = dot_n(data.u));

[error_u, error_p, error_phat, h_list, N_list] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    N    = 2^i;
    mesh = StructuredSimplicialGrid(N, N);
    mesh = apply_tag(mesh, "Neumann", @(X) isapprox(abs(X(1, :) - 0.5), 0.5));

    solution = solve(problem, mesh, Order= k);
    [error, error_hat] = solution.compute_l2_error({data.u, data.p}, {data.p});

    N_list(i)     = N;
    h_list(i)     = sqrt(2) / N;
    error_u(i)    = error(1);
    error_p(i)    = error(2);
    error_phat(i) = error_hat(1);
end

r_u    = [NaN; compute_eoc(error_u,    h_list)];
r_p    = [NaN; compute_eoc(error_p,    h_list)];
r_phat = [NaN; compute_eoc(error_phat, h_list)];

hdr = "   N  |     h    |   e(u_h)  | r(u_h) |   e(p_h)  | r(p_h) |  e(p_hat) | r(p_hat)";
fprintf("\n%s\n%s\n", hdr, repmat('-', 1, strlength(hdr)));
for i = 1:max_iterations
    if i == 1
        fprintf(" %4d | %8.2e | %9.3e |   -    | %9.3e |   -    | %9.3e |   -\n", ...
            N_list(i), h_list(i), error_u(i), error_p(i), error_phat(i));
    else
        fprintf(" %4d | %8.2e | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f\n", ...
            N_list(i), h_list(i), error_u(i), r_u(i), error_p(i), r_p(i), error_phat(i), r_phat(i));
    end
end
end
