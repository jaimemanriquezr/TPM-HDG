%-------------------------------------------------------------------------
%Description: The Cockburn–Shen stabilization of HyperelasticGradientProblem
%             (stabilization = "cockburn_shen": τ_KF = C_KF,o/h_F on
%             T_h^(-,nb), τ = 1 elsewhere; Cockburn and Shen, Results in
%             Applied Mathematics 1 (2019), Thm. 2.4, Sec. 2.5.5, Sec. 5),
%             uniform degree k (their degree raise on T_h^(-,nb) is not
%             implemented), τ lagged.
%             (A) Their expansion test (Sec. 5.2.1, eq. 5.2): neo-Hookean,
%                 E = 1, ν = 0.25, u = (0.1 (1−Y) Y cos(πX), Y), Dirichlet
%                 data, N = 2, 4, 8, 16, 4 load steps. A(F_h) has no negative
%                 eigenvalue, so τ = 1 everywhere, and F_h, u_h, û_h converge
%                 at k + 1 (their Table 2).
%             (B) Their compression test (Sec. 5.2.2, eq. 5.3), k = 1,
%                 N = 8, 16, 32, Newton started from the projection of the
%                 exact solution (from rest our Newton, without their line
%                 search, does not converge): u_h at 2, F_h at 1 (their
%                 Table 5).
%             (C) The manufactured problems of
%                 test_hyperelastic_gradient_convergence (u = a (sin(πx) y,
%                 x cos(πy)), E = 1, ν = 0.3): Saint Venant–Kirchhoff and
%                 neo-Hookean at a = 0.05, neo-Hookean at a = 0.1, k = 1, 2,
%                 N = 4, 8, 16, from rest. Prints the EOC of F_h, P(F_h), u_h,
%                 û_h, the Newton steps, the fraction of elements in
%                 T_h^(-,nb), and the range of τ h on them against
%                 40 μ h/h_F.
%             Checks: (A) EOC ≥ k + 0.8 in all three fields and no element
%             marked; (B) EOC u ≥ 1.8, F ≥ 0.8; (C) Newton converged, EOC
%             u, û ≥ k + 0.8 and F, P ≥ k − 0.2.
%-------------------------------------------------------------------------
n_failed = 0;
warning("off", "HyperelasticGradientProblem:newtonNotConverged");

% (A) expansion
E = 1; nu = 0.25;
data = HyperelasticGradientProblem.get_manufactured_solution(["0.1*(1-y)*y*cos(pi*x)"; "y"], E, nu, 0, 0, 0, "neo_hookean");
Ns = [2 4 8 16];
fprintf("(A) expansion test, neo-Hookean, tau = cockburn_shen\n");
for k = 1:2
    [e, marked] = deal(zeros(numel(Ns), 3), zeros(numel(Ns), 1));
    for i = 1:numel(Ns)
        mesh = StructuredSimplicialGrid(Ns(i), Ns(i));
        pde = HyperelasticGradientProblem(k, material = "neo_hookean", young_modulus = E, poisson_ratio = nu, ...
            source_term = data.f, datum_dirichlet = data.u, stabilization = "cockburn_shen", load_steps = 4);
        [solution, ~, newton] = solve(pde, mesh, k);
        [ev, eh] = solution.compute_l2_error({data.gradient, data.u}, {data.u});
        e(i, :) = [ev(:).', eh];
        marked(i) = mean(newton.marked);
    end
    r = compute_eoc(e, 1./Ns(:));
    fprintf("  k = %d: e(F_h) = %s; EOC F, u, û (last) = %s; marked %s\n", k, ...
        sprintf("%.2e ", e(:, 1)), sprintf("%.2f ", r(end, :)), sprintf("%.2f ", marked));
    n_failed = n_failed + check(all(r(end, :) >= k + 0.8) && all(marked == 0), ...
        sprintf("(A) k = %d: tau = 1 and order k + 1", k));
end

% (B) compression, k = 1, from the projected exact solution
data = HyperelasticGradientProblem.get_manufactured_solution( ...
    ["-0.04375*(1-y)*y*cos(pi*x)"; "0.15*y-0.12*y^2*(1-cos(2*pi*x))"], E, nu, 0, 0, 0, "neo_hookean");
Ns = [8 16 32]; k = 1;
e = zeros(numel(Ns), 3); marked = zeros(numel(Ns), 1);
for i = 1:numel(Ns)
    mesh = StructuredSimplicialGrid(Ns(i), Ns(i));
    G = data.gradient;
    iterate = struct("gradient", [L2proj2D(mesh, G{1, 1}, k); L2proj2D(mesh, G{2, 1}, k);
        L2proj2D(mesh, G{1, 2}, k); L2proj2D(mesh, G{2, 2}, k)], ...
        "displacement", L2proj2D(mesh, data.u, k), "displacement_skeleton", L2projSkeleton1d(mesh, data.u, k));
    pde = HyperelasticGradientProblem(k, material = "neo_hookean", young_modulus = E, poisson_ratio = nu, ...
        source_term = data.f, datum_dirichlet = data.u, stabilization = "cockburn_shen", iterate = iterate);
    [solution, ~, newton] = solve(pde, mesh, k);
    [ev, eh] = solution.compute_l2_error({data.gradient, data.u}, {data.u});
    e(i, :) = [ev(:).', eh];
    marked(i) = mean(newton.marked);
end
r = compute_eoc(e, 1./Ns(:));
fprintf("(B) compression test, k = 1: e(u_h) = %s, e(F_h) = %s; EOC u = %s, F = %s; marked %s\n", ...
    sprintf("%.2e ", e(:, 2)), sprintf("%.2e ", e(:, 1)), sprintf("%.2f ", r(:, 2)), ...
    sprintf("%.2f ", r(:, 1)), sprintf("%.2f ", marked));
n_failed = n_failed + check(r(end, 2) >= 1.8 && r(end, 1) >= 0.8, "(B) orders 2 (u) and 1 (F)");

% (C) the manufactured problems of this branch
E = 1; nu = 0.3;
[~, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
cases = {"saint_venant_kirchhoff", 0.05; "neo_hookean", 0.05; "neo_hookean", 0.1};
Ns = [4 8 16];
fprintf("(C) a (sin(pi x) y, x cos(pi y)), tau = cockburn_shen (40 mu h/h_F = %.1f to %.1f)\n", ...
    40*mu/sqrt(2), 40*mu);
for c = 1:size(cases, 1)
    [material, a] = cases{c, :};
    data = HyperelasticGradientProblem.get_manufactured_solution( ...
        [a + "*sin(pi*x)*y"; a + "*x*cos(pi*y)"], E, nu, 0, 0, 0, material);
    for k = 1:2
        [e, its, converged, marked, tau_h] = deal(zeros(numel(Ns), 4), zeros(numel(Ns), 1), ...
            false(numel(Ns), 1), zeros(numel(Ns), 1), zeros(numel(Ns), 2));
        for i = 1:numel(Ns)
            mesh = StructuredSimplicialGrid(Ns(i), Ns(i));
            pde = HyperelasticGradientProblem(k, material = material, young_modulus = E, poisson_ratio = nu, ...
                source_term = data.f, datum_dirichlet = data.u, stabilization = "cockburn_shen");
            [solution, ~, newton] = solve(pde, mesh, k);
            [ev, eh] = solution.compute_l2_error({data.gradient, data.u}, {data.u});
            e(i, :) = [ev(1), piola_error(solution, data, pde, mesh, k), ev(2), eh];
            its(i) = newton.iterations;
            converged(i) = newton.converged;
            marked(i) = mean(newton.marked);
            t = newton.tau(:, newton.marked)/Ns(i);
            if isempty(t), t = 0; end
            tau_h(i, :) = [min(t, [], 'all'), max(t, [], 'all')];
        end
        r = compute_eoc(e, 1./Ns(:));
        fprintf("  %s a = %g k = %d: EOC F, P, u, û = %s; steps %s; marked %s; tau h on marked in [%.0f, %.0f]\n", ...
            material, a, k, sprintf("%.2f ", r(end, :)), sprintf("%d ", its), sprintf("%.2f ", marked), ...
            min(tau_h(:, 1)), max(tau_h(:, 2)));
        tag = sprintf("(C) %s a = %g k = %d: ", material, a, k);
        n_failed = n_failed + check(all(converged), tag + "Newton converged");
        n_failed = n_failed + check(all(r(end, 3:4) >= k + 0.8) && all(r(end, 1:2) >= k - 0.2), tag + "orders");
    end
end
fprintf("\n%d check(s) failed.\n", n_failed);

function err = piola_error(solution, data, pde, mesh, k)
% L² error of P(I + H_h) against the exact P, by quadrature
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(pde.Data.young_modulus, pde.Data.poisson_ratio);
[xy, rs, w] = quadpoints_elements_2d(mesh, k + 3, output="cell");
p = dubiner_basis(k, rs{:});
H = solution.Unknowns.field(1).dof;
[n_q, d2] = size(p);
N = size(H, 2);
G = pagemtimes(p, reshape(H, d2, 4, N));
P = HyperelasticGradientProblem.first_piola(reshape(permute(G, [1 3 2]), [], 4), lambda, mu, pde.Data.material);
P = reshape(P, n_q, N, 4);
P_exact = reshape(cat(3, data.piola{1, 1}(xy{:}), data.piola{2, 1}(xy{:}), ...
    data.piola{1, 2}(xy{:}), data.piola{2, 2}(xy{:})), n_q, N, 4);
err = sqrt(sum(reshape(mesh.volume, 1, []) .* (w.' * sum((P - P_exact).^2, 3))));
end
