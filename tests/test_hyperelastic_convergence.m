%-------------------------------------------------------------------------
%Description: Spatial convergence of HyperelasticProblem (steady).
%
%             Manufactured displacement
%               u = a (sin(πx) y, x cos(πy)),  a = 0.05,
%             on the unit square (E = 1, ν = 0.3, Dirichlet data on the
%             whole boundary): |∇u| reaches 0.16 and the stress remainder
%             N(∇u) is about 8% of the linear stress. Structured meshes
%             N = 4, 8, 16. Errors against the exact first Piola stress
%             σ = P (field 1), u, γ = skew ∇u and û.
%             Two stabilizations (mesh.tau):
%               "h"  τ = 4/h (about 10μ/h): the robust choice; u_h, û_h
%                    converge at k + 1, σ_h, γ_h at k (the O(1/h)
%                    stabilization costs the flux one order, as in the
%                    linear HDG);
%               "1"  τ = 1 (the linear class's default): all four fields
%                    at k + 1, shown for Saint Venant–Kirchhoff, k = 2 (at
%                    k = 1, or on finer meshes, τ = O(1) lets the discrete
%                    tangent lose its control of rotations).
%             Checks (EOC of the last two meshes):
%               (1) τ = 4/h, Saint Venant–Kirchhoff k = 1, 2 and
%                   neo-Hookean k = 2: EOC u_h, û_h ≥ k + 0.8 and
%                   EOC σ_h, γ_h ≥ k − 0.2;
%               (2) τ = 1, Saint Venant–Kirchhoff k = 2: EOC of all four
%                   fields ≥ k + 0.8;
%               (3) the linear solve (LinearElasticityProblem with the same
%                   source, data and mesh) stalls: its u_h error on the
%                   finest mesh exceeds 10x the hyperelastic one;
%               (4) Newton converged on every mesh.
%-------------------------------------------------------------------------
n_failed = 0;
E = 1; nu = 0.3;
u_exact = ["0.05*sin(pi*x)*y"; "0.05*x*cos(pi*y)"];
Ns = [4 8 16];
cases = struct("material", {"saint_venant_kirchhoff", "neo_hookean", "saint_venant_kirchhoff", "saint_venant_kirchhoff"}, ...
    "k", {1, 2, 2, 2}, "tau", {"h", "h", "h", "1"});

for c = cases
    data = HyperelasticProblem.get_manufactured_solution(u_exact, E, nu, 0, 0, 0, c.material);
    pde = HyperelasticProblem(c.k, material = c.material, young_modulus = E, poisson_ratio = nu, ...
        source_term = data.f, datum_dirichlet = data.u);
    lin = LinearElasticityProblem(c.k, young_modulus = E, poisson_ratio = nu, ...
        source_term = data.f, datum_dirichlet = data.u);
    n = numel(Ns);
    [e, e_lin, iterations, converged] = deal(zeros(n, 4), zeros(n, 1), zeros(n, 1), false(n, 1));
    for i = 1:n
        mesh = StructuredSimplicialGrid(Ns(i), Ns(i));
        if c.tau == "h", mesh = mesh.add_stabilization(4*Ns(i)); end
        [solution, ~, newton] = solve(pde, mesh, c.k);
        iterations(i) = newton.iterations;
        converged(i) = newton.converged;
        [e_vol, e_hat] = solution.compute_l2_error({data.sigma, data.u, data.gamma}, {data.u});
        e(i, :) = [e_vol(:).', e_hat];
        e_vol_lin = solve(lin, mesh, c.k).compute_l2_error({data.sigma, data.u, data.gamma}, {data.u});
        e_lin(i) = e_vol_lin(2);
    end
    h = 1./Ns(:);
    r = [nan(1, 4); compute_eoc(e, h)];
    r_lin = [nan; compute_eoc(e_lin, h)];

    fprintf("\n%s, k = %d, tau = %s\n", c.material, c.k, c.tau);
    fprintf("   N | it |  e(σ_h)    EOC |  e(u_h)    EOC |  e(γ_h)    EOC |  e(û_h)    EOC | linear e(u_h)  EOC\n");
    for i = 1:n
        fprintf(" %3d | %2d | %.3e %5.2f | %.3e %5.2f | %.3e %5.2f | %.3e %5.2f |   %.3e %5.2f\n", ...
            Ns(i), iterations(i), e(i, 1), r(i, 1), e(i, 2), r(i, 2), e(i, 3), r(i, 3), ...
            e(i, 4), r(i, 4), e_lin(i), r_lin(i));
    end
    tag = sprintf("%s k = %d tau = %s: ", c.material, c.k, c.tau);
    names = ["sigma", "u", "gamma", "u_hat"];
    target = (c.k + 0.8)*[1 1 1 1];
    if c.tau == "h", target([1 3]) = c.k - 0.2; end
    for j = 1:4
        n_failed = n_failed + check(r(end, j) >= target(j), ...
            tag + "EOC " + names(j) + " >= " + target(j));
    end
    n_failed = n_failed + check(e_lin(end) > 10*e(end, 2), tag + "the linear solve stalls");
    n_failed = n_failed + check(all(converged), tag + "Newton converged");
end
fprintf("\n%d check(s) failed.\n", n_failed);
