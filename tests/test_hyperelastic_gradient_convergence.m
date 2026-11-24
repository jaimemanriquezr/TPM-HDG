%-------------------------------------------------------------------------
%Description: Spatial convergence of HyperelasticGradientProblem (steady).
%
%             Manufactured displacement u = a (sin(πx) y, x cos(πy)) on
%             the unit square (E = 1, ν = 0.3, Dirichlet data on the whole
%             boundary): Saint Venant–Kirchhoff at a = 0.05, neo-Hookean at
%             a = 0.05 and 0.1 (Saint Venant–Kirchhoff loses strong
%             ellipticity on this solution from a ≈ 0.075, so a = 0.1 is
%             not a valid test for it); k = 1, 2, structured meshes
%             N = 4, 8, 16, default stabilization τ = 40 μ / h_F, Newton
%             from rest.
%             Errors (L²) of F_h = I + H_h (field 1), P(F_h) (quadrature),
%             u_h and û_h; EOC of the last two meshes.
%             Checks:
%               (1) Newton converged on every mesh;
%               (2) EOC u_h, û_h ≥ k + 0.8;
%               (3) EOC F_h, P(F_h) ≥ k − 0.2 (with τ = O(1/h) the
%                   gradient converges at order k, as in the linear HDG of
%                   the same form).
%             Set the variable cases (rows {material, a}) or load_steps
%             before running to change them.
%-------------------------------------------------------------------------
n_failed = 0;
if ~exist("cases", "var")
    cases = {"saint_venant_kirchhoff", 0.05; "neo_hookean", 0.05; "neo_hookean", 0.1};
end
if ~exist("load_steps", "var"), load_steps = 1; end
E = 1; nu = 0.3; Ns = [4 8 16];
for c = 1:size(cases, 1)
    [material, a] = cases{c, :};
    u_exact = [a + "*sin(pi*x)*y"; a + "*x*cos(pi*y)"];
    data = HyperelasticGradientProblem.get_manufactured_solution(u_exact, E, nu, 0, 0, 0, material);
    for k = 1:2
        pde = HyperelasticGradientProblem(k, material = material, young_modulus = E, ...
            poisson_ratio = nu, source_term = data.f, datum_dirichlet = data.u, load_steps = load_steps);
        [e, its, converged] = deal(zeros(numel(Ns), 4), zeros(numel(Ns), 1), false(numel(Ns), 1));
        for i = 1:numel(Ns)
            mesh = StructuredSimplicialGrid(Ns(i), Ns(i));
            [solution, ~, newton] = solve(pde, mesh, k);
            [e_vol, e_hat] = solution.compute_l2_error({data.gradient, data.u}, {data.u});
            e(i, :) = [e_vol(1), piola_error(solution, data, pde, mesh, k), e_vol(2), e_hat];
            its(i) = newton.iterations;
            converged(i) = newton.converged;
        end
        r = [nan(1, 4); compute_eoc(e, 1./Ns(:))];
        fprintf("\n%s, a = %g, k = %d\n   N | it |  e(F_h)    EOC |  e(P_h)    EOC |  e(u_h)    EOC |  e(û_h)    EOC\n", ...
            material, a, k);
        for i = 1:numel(Ns)
            fprintf(" %3d | %2d | %.3e %5.2f | %.3e %5.2f | %.3e %5.2f | %.3e %5.2f\n", Ns(i), its(i), ...
                e(i, 1), r(i, 1), e(i, 2), r(i, 2), e(i, 3), r(i, 3), e(i, 4), r(i, 4));
        end
        tag = sprintf("%s a = %g k = %d: ", material, a, k);
        n_failed = n_failed + check(all(converged), tag + "Newton converged");
        n_failed = n_failed + check(all(r(end, 3:4) >= k + 0.8), tag + "EOC u, u_hat >= k + 0.8");
        n_failed = n_failed + check(all(r(end, 1:2) >= k - 0.2), tag + "EOC F, P >= k - 0.2");
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
