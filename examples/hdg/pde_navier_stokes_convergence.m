function pde_navier_stokes_convergence(k, max_iterations, formulation, stream_function, pressure, viscosity)
% Manufactured steady Navier–Stokes by Picard on the Oseen solver
% (Solano & Vargas 2022, §6.4), on either fluid formulation:
%   "gradient": OseenProblem        (L = grad u, p)        errors in L, u, p
%   "stress":   NavierStokesProblem (PEERS sigma, gamma)   errors in sigma, u, gamma
% Structured meshes of the unit square; prints errors, EOC and Picard
% iteration counts. Expected: EOC ~ k+1, and ~ k+2 for the divergence-free u*.
arguments
    k               = 2
    max_iterations  = 5
    formulation     = "gradient"
    stream_function = "sin(pi*x)^2*sin(pi*y)^2"
    pressure        = "cos(pi*x)*sin(pi*y)"
    viscosity       = 1e-1
end

switch formulation
    case "gradient"
        data    = OseenProblem.get_manufactured_solution(stream_function, pressure, viscosity, []);
        problem = OseenProblem(viscosity = viscosity, source_term = data.f, datum_dirichlet = data.u);
        exact   = {data.L, data.u, data.p};   names = ["L_h", "u_h", "p_h"];
    case "stress"
        data    = NavierStokesProblem.get_manufactured_solution(stream_function, pressure, viscosity, []);
        problem = NavierStokesProblem(viscosity = viscosity, source_term = data.f, datum_dirichlet = data.u);
        exact   = {data.sigma, data.u, data.gamma};   names = ["\sigma_h", "u_h", "\gamma_h"];
end

[errors, h_list, N_list, iterations] = deal(zeros(max_iterations, 4), zeros(max_iterations, 1), zeros(max_iterations, 1), zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    N    = 2^(i + 1);
    mesh = StructuredSimplicialGrid(N, N);
    [solution, history, u_star] = solve_navier_stokes(problem, mesh, k, verbose = true);
    error = solution.compute_l2_error(exact);
    error(4) = l2_error_p1(mesh, k, u_star, data.u);   % divergence-free u*, expected k+2

    N_list(i)     = N;
    h_list(i)     = sqrt(2) / N;
    errors(i, :)  = error(1:4);
    iterations(i) = numel(history);
end

rates = [NaN(1, 4); log(errors(1:end-1, :)./errors(2:end, :)) ./ log(h_list(1:end-1)./h_list(2:end))];
hdr = sprintf("   N  |     h    | Picard |  e(%s)  |  r  |  e(%s)  |  r  |  e(%s)  |  r  |  e(u*)  |  r", names);
fprintf("\n%s\n%s\n", hdr, repmat('-', 1, strlength(hdr)));
for i = 1:max_iterations
    fprintf(" %4d | %8.2e |  %3d   | %9.3e | %4.2f | %9.3e | %4.2f | %9.3e | %4.2f | %9.3e | %4.2f\n", ...
        N_list(i), h_list(i), iterations(i), errors(i, 1), rates(i, 1), errors(i, 2), rates(i, 2), errors(i, 3), rates(i, 3), errors(i, 4), rates(i, 4));
end
end

function e = l2_error_p1(mesh, k, u_star, u_exact)
% L2 error of a P_{k+1} Dubiner vector field against handles
[xy, rs, w] = quadpoints_elements_2d(mesh, k + 2, output = "cell");
P1 = dubiner_basis(k + 1, rs{:});
e = 0;
for i = 1:2
    exact = u_exact{i}(xy{:}) + 0*xy{1};
    e = e + sum(mesh.volume .* sum(w .* (P1*u_star{i} - exact).^2, 1));
end
e = sqrt(e);
end
