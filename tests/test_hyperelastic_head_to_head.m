%-------------------------------------------------------------------------
%Description: HyperelasticGradientProblem (gradient-based HDG) against
%             HyperelasticProblem (PEERS-based HDG) on the same data.
%
%             u = a (sin(πx) y, x cos(πy)), E = 1, ν = 0.3, Dirichlet data,
%             N = 8, 16, k = 1, 2, Newton from rest. Methods:
%               PEERS  τ = 1      HyperelasticProblem, mesh default
%               PEERS  τ = 4/h    HyperelasticProblem
%               grad   τ = 10μ    HyperelasticGradientProblem, "constant", ℓ = 1
%               grad   τ = 40μ/h  HyperelasticGradientProblem default
%             Cases: Saint Venant–Kirchhoff a = 0.05, neo-Hookean a = 0.05
%             and 0.1 (set the variable cases to choose). Prints per
%             method the Newton steps ("!" not converged, "x" inverted
%             element), the errors of u_h, û_h and of the first Piola
%             stress (PEERS: σ_h; gradient: P(I + H_h)) on N = 16 and
%             their EOC from N = 8.
%             Checks (the default gradient method):
%               (1) Newton converged in at most 6 steps everywhere;
%               (2) EOC u_h, û_h ≥ k + 0.8 and EOC P ≥ k − 0.2.
%-------------------------------------------------------------------------
n_failed = 0;
if ~exist("cases", "var")
    cases = {"saint_venant_kirchhoff", 0.05; "neo_hookean", 0.05; "neo_hookean", 0.1};
end
E = 1; nu = 0.3; Ns = [8 16];
methods = ["PEERS tau=1", "PEERS tau=4/h", "grad tau=10mu", "grad tau=40mu/h"];
warning("off", "HyperelasticProblem:newtonNotConverged");
warning("off", "HyperelasticGradientProblem:newtonNotConverged");
for c = 1:size(cases, 1)
    [material, a] = cases{c, :};
    u_exact = [a + "*sin(pi*x)*y"; a + "*x*cos(pi*y)"];
    data_g = HyperelasticGradientProblem.get_manufactured_solution(u_exact, E, nu, 0, 0, 0, material);
    data_p = HyperelasticProblem.get_manufactured_solution(u_exact, E, nu, 0, 0, 0, material);
    for k = 1:2
        fprintf("\n%s, a = %g, k = %d (N = 16; EOC from N = 8)\n", material, a, k);
        fprintf("  %-16s | steps  |  e(u_h)    EOC  |  e(û_h)    EOC  |  e(P_h)    EOC\n", "method");
        for m = 1:numel(methods)
            [e, steps] = deal(nan(numel(Ns), 3), strings(1, numel(Ns)));
            for i = 1:numel(Ns)
                mesh = StructuredSimplicialGrid(Ns(i), Ns(i));
                try
                    switch m
                        case {1, 2}
                            if m == 2, mesh = mesh.add_stabilization(4*Ns(i)); end
                            pde = HyperelasticProblem(k, material = material, young_modulus = E, ...
                                poisson_ratio = nu, source_term = data_p.f, datum_dirichlet = data_p.u);
                            [solution, ~, newton] = solve(pde, mesh, k);
                            [ev, eh] = solution.compute_l2_error({data_p.sigma, data_p.u, data_p.gamma}, {data_p.u});
                            e(i, :) = [ev(2), eh, ev(1)];
                        otherwise
                            options = {};
                            if m == 3, options = {"stabilization", "constant", "stabilization_constant", 10}; end
                            pde = HyperelasticGradientProblem(k, "material", material, "young_modulus", E, ...
                                "poisson_ratio", nu, "source_term", data_g.f, "datum_dirichlet", data_g.u, options{:});
                            [solution, ~, newton] = solve(pde, mesh, k);
                            [ev, eh] = solution.compute_l2_error({data_g.gradient, data_g.u}, {data_g.u});
                            e(i, :) = [ev(2), eh, piola_error(solution, data_g, pde, mesh, k)];
                    end
                    steps(i) = string(newton.iterations) + string(repmat('!', 1, ~newton.converged));
                catch err
                    if ~startsWith(err.identifier, "HyperelasticProblem:inverted"), rethrow(err); end
                    steps(i) = "x";
                end
            end
            r = log(e(1, :)./e(2, :))/log(2);
            fprintf("  %-16s | %-6s | %.2e %5.2f  | %.2e %5.2f  | %.2e %5.2f\n", methods(m), ...
                strjoin(steps, ","), e(2, 1), r(1), e(2, 2), r(2), e(2, 3), r(3));
            if m == 4
                tag = sprintf("%s a = %g k = %d: ", material, a, k);
                ok = all(~contains(steps, ["!", "x"])) && all(double(steps) <= 6);
                n_failed = n_failed + check(ok, tag + "gradient default converged in <= 6 steps");
                n_failed = n_failed + check(all(r(1:2) >= k + 0.8) && r(3) >= k - 0.2, ...
                    tag + "gradient default EOC");
            end
        end
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
