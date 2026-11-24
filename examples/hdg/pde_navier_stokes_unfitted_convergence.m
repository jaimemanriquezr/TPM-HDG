function pde_navier_stokes_unfitted_convergence(k, max_iterations, formulation, viscosity)
% Unfitted steady Navier–Stokes on the disk of radius 0.75 immersed in a
% background grid of [-1,1]^2 (Solano & Vargas 2022, Example 4 geometry,
% Example 3 construction: the boundary Gamma_h lies at distance O(h) from
% the circle). The Dirichlet datum is transferred radially along the
% transfer paths (transfer_mapping), Navier–Stokes by Picard on the Oseen
% solver with beta = u*. Expected: k+1 for L/sigma, u, p and about k+2 for
% the divergence-free u* (their Tables 3-4).
%   "gradient": OseenProblem (L = grad u, p);  "stress": NavierStokesProblem
arguments
    k               = 2
    max_iterations  = 4
    formulation     = "gradient"
    viscosity       = 1
end
R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial = @(X) R * X ./ vecnorm(X);
% the paper's Example 4 solution (here p = sin(x^2 + y^2), mean removed)
U = "-cos(x)*sin(y)";            % u = (-dU/dy, dU/dx) = (cos x cos y, sin x sin y)
P = "sin(x^2 + y^2)";

switch formulation
    case "gradient"
        data    = OseenProblem.get_manufactured_solution(U, P, viscosity, []);
        problem = OseenProblem(viscosity = viscosity, source_term = data.f, datum_dirichlet = data.u, transfer_mapping = radial);
        exact   = {data.L, data.u, data.p};   names = ["L_h", "u_h", "p_h"];
    case "stress"
        data    = NavierStokesProblem.get_manufactured_solution(U, P, viscosity, []);
        problem = NavierStokesProblem(viscosity = viscosity, source_term = data.f, datum_dirichlet = data.u, transfer_mapping = radial);
        exact   = {data.sigma, data.u, data.gamma};   names = ["\sigma_h", "u_h", "\gamma_h"];
end

[errors, h_list, N_list, iterations] = deal(zeros(max_iterations, 4), zeros(max_iterations, 1), zeros(max_iterations, 1), zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    N    = 8 * 2^(i - 1);
    mesh = immerse_mesh(StructuredSimplicialGrid(N, N, -1, 1, -1, 1), phi);
    [solution, history, u_star] = solve_navier_stokes(problem, mesh, k, verbose = true);
    % the discrete pressure has zero mean on D_h, not on the unit square:
    % shift the exact p (and sigma = 2 mu eps - p I) by that constant
    if formulation == "gradient"
        exact{3} = zero_mean_on(mesh, k, data.p);
    else
        [~, c] = zero_mean_on(mesh, k, @(x, y) -(data.sigma{1, 1}(x, y) + data.sigma{2, 2}(x, y))/2);
        exact{1} = data.sigma;
        exact{1}{1, 1} = @(x, y) data.sigma{1, 1}(x, y) + c;
        exact{1}{2, 2} = @(x, y) data.sigma{2, 2}(x, y) + c;
    end
    error = solution.compute_l2_error(exact);
    error(4) = l2_error_p1(mesh, k, u_star, data.u);

    N_list(i)     = N;
    h_list(i)     = 2 / N;
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

function [p0, c] = zero_mean_on(mesh, k, p)
[xy, ~, w] = quadpoints_elements_2d(mesh, k + 2, output = "cell");
c = sum(mesh.volume .* sum(w .* p(xy{:}), 1)) / sum(mesh.volume);
p0 = @(x, y) p(x, y) - c;
end

function e = l2_error_p1(mesh, k, u_star, u_exact)
[xy, rs, w] = quadpoints_elements_2d(mesh, k + 2, output = "cell");
P1 = dubiner_basis(k + 1, rs{:});
e = 0;
for i = 1:2
    exact = u_exact{i}(xy{:}) + 0*xy{1};
    e = e + sum(mesh.volume .* sum(w .* (P1*u_star{i} - exact).^2, 1));
end
e = sqrt(e);
end
