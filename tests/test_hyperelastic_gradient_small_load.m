%-------------------------------------------------------------------------
%Description: Small-load limit of HyperelasticGradientProblem.
%
%             Load s (f_s, g_s) = s (f, g), (f, g) the linear manufactured
%             data of the displacement (sin(πx) y, x cos(πy)) (E = 1,
%             ν = 0.3), N = 4, k = 2, default stabilization, both
%             materials, s = 0.04, 0.02, 0.01, 0.005, against the linear
%             HDG elasticity of the same form (material = "linear", the
%             same spaces, flux and τ).
%             Checks:
%               (1) at the zero iterate the hyperelastic operators equal
%                   the linear ones (relative difference ≤ 1e-12);
%               (2) the difference of all local dofs (H, u) and of û has
%                   order 2 in s (EOC ≥ 1.9);
%               (3) the relative difference is below 5 s.
%-------------------------------------------------------------------------
n_failed = 0;
E = 1; nu = 0.3; k = 2;
mesh = StructuredSimplicialGrid(4, 4);
data = HyperelasticGradientProblem.get_manufactured_solution(["sin(pi*x)*y"; "x*cos(pi*y)"], E, nu, 0, 0, 0, "linear");
scaled = @(h, s) {@(x, y) s*h{1}(x, y); @(x, y) s*h{2}(x, y)};
loads = [0.04 0.02 0.01 0.005];
make = @(material, s) HyperelasticGradientProblem(k, material = material, young_modulus = E, ...
    poisson_ratio = nu, source_term = scaled(data.f, s), datum_dirichlet = scaled(data.u, s));

tau_mesh = mesh.add_stabilization(1);
tau_mesh.tau = 40*0.5*E/(1 + nu)./mesh.area(mesh.boundaries);
ops = get_matrices(make("saint_venant_kirchhoff", 1), tau_mesh, k);
ops_lin = get_matrices(make("linear", 1), tau_mesh, k);
d = max(cellfun(@(A, B) max(abs(A - B), [], 'all')/max(abs(B), [], 'all'), ops, ops_lin));
fprintf("zero iterate: max relative difference of the local operators %.2e\n", d);
n_failed = n_failed + check(d <= 1e-12, "zero iterate: linear operators");

for material = ["saint_venant_kirchhoff", "neo_hookean"]
    [d_local, d_trace, scale] = deal(zeros(size(loads)));
    for i = 1:numel(loads)
        [solution, ~, newton] = solve(make(material, loads(i)), mesh, k);
        solution_lin = solve(make("linear", loads(i)), mesh, k);
        X = [solution.Unknowns.field(1).dof; cat(1, solution.Unknowns.field(2).dof)];
        X_lin = [solution_lin.Unknowns.field(1).dof; cat(1, solution_lin.Unknowns.field(2).dof)];
        d_local(i) = max(abs(X - X_lin), [], 'all');
        d_trace(i) = max(abs(cat(1, solution.Unknowns.hybrid.dof) - cat(1, solution_lin.Unknowns.hybrid.dof)), [], 'all');
        scale(i) = max(abs(X_lin), [], 'all');
        n_failed = n_failed + check(newton.converged, material + ": Newton converged");
    end
    eoc = @(v) [nan, log(v(1:end-1)./v(2:end))./log(loads(1:end-1)./loads(2:end))];
    r_local = eoc(d_local); r_trace = eoc(d_trace);
    fprintf("\n%s\n      s    | max|X - X_lin|  EOC  | max|û - û_lin|  EOC  | relative\n", material);
    for i = 1:numel(loads)
        fprintf("  %.4f  |   %.3e   %5.2f |   %.3e   %5.2f | %.2e\n", loads(i), d_local(i), ...
            r_local(i), d_trace(i), r_trace(i), d_local(i)/scale(i));
    end
    n_failed = n_failed + check(all(r_local(2:end) >= 1.9) && all(r_trace(2:end) >= 1.9), ...
        material + ": difference of order load^2");
    n_failed = n_failed + check(all(d_local./scale <= 5*loads), material + ": relative difference below 5 load");
end
fprintf("\n%d check(s) failed.\n", n_failed);
