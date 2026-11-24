%-------------------------------------------------------------------------
%Description: Small-load limit of HyperelasticProblem.
%
%             Load s (f_s, g_s) = s (f, g), with (f, g) the source and
%             Dirichlet data of the linear manufactured displacement
%             (sin(πx) y, x cos(πy)) (LinearElasticityProblem
%             .get_manufactured_solution, E = 1, ν = 0.3), on a structured
%             mesh N = 4, k = 2, τ = 1, both materials, s = 0.04, 0.02,
%             0.01, 0.005. The hyperelastic solution X_s (all local dofs
%             [σ; u; γ] and û) and the linear one X_s^lin = s X_1^lin
%             differ by the stress remainder, quadratic at zero:
%             |X_s − X_s^lin| = O(s²).
%             Checks:
%               (1) with no iterate, get_matrices of the hyperelastic class
%                   returns exactly the linear class's operators;
%               (2) the difference d(s) = max |X_s − X_s^lin| has order 2
%                   in s (EOC ≥ 1.9 between consecutive s) for the local
%                   dofs and the trace;
%               (3) the relative difference d(s)/max|X_s^lin| is below
%                   5 s (|∇u_s| reaches πs; the linear class is the
%                   small-load limit).
%-------------------------------------------------------------------------
n_failed = 0;
E = 1; nu = 0.3; k = 2;
mesh = StructuredSimplicialGrid(4, 4);
data = LinearElasticityProblem.get_manufactured_solution(["sin(pi*x)*y"; "x*cos(pi*y)"], E, nu);
scaled = @(h, s) {@(x, y) s*h{1}(x, y); @(x, y) s*h{2}(x, y)};
loads = [0.04 0.02 0.01 0.005];

% (1) zero iterate: the linear operators
pde = HyperelasticProblem(k, young_modulus = E, poisson_ratio = nu, source_term = data.f, datum_dirichlet = data.u);
lin = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, source_term = data.f, datum_dirichlet = data.u);
[ops, lhs, rhs] = get_matrices(pde, mesh, k);
[ops_lin, lhs_lin, rhs_lin] = get_matrices(lin, mesh, k);
n_failed = n_failed + check(isequal(ops, ops_lin) && isequal(lhs, lhs_lin) && isequal(rhs, rhs_lin), ...
    "zero iterate: get_matrices is the linear class's");

for material = ["saint_venant_kirchhoff", "neo_hookean"]
    [d_local, d_trace, scale] = deal(zeros(size(loads)));
    for i = 1:numel(loads)
        s = loads(i);
        pde = HyperelasticProblem(k, material = material, young_modulus = E, poisson_ratio = nu, ...
            source_term = scaled(data.f, s), datum_dirichlet = scaled(data.u, s));
        lin = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, ...
            source_term = scaled(data.f, s), datum_dirichlet = scaled(data.u, s));
        [solution, ~, newton] = solve(pde, mesh, k);
        solution_lin = solve(lin, mesh, k);
        X = local_dofs(solution);
        X_lin = local_dofs(solution_lin);
        x_hat = cat(1, solution.Unknowns.hybrid.dof);
        x_hat_lin = cat(1, solution_lin.Unknowns.hybrid.dof);
        d_local(i) = max(abs(X - X_lin), [], 'all');
        d_trace(i) = max(abs(x_hat - x_hat_lin), [], 'all');
        scale(i) = max(abs(X_lin), [], 'all');
        n_failed = n_failed + check(newton.converged, material + ": Newton converged");
    end
    eoc = @(d) log(d(1:end-1)./d(2:end))./log(loads(1:end-1)./loads(2:end));
    fprintf("\n%s\n      s    | max|X - X_lin|  EOC  | max|û - û_lin|  EOC  | relative\n", material);
    r_local = [nan, eoc(d_local)];
    r_trace = [nan, eoc(d_trace)];
    for i = 1:numel(loads)
        fprintf("  %.4f  |   %.3e   %5.2f |   %.3e   %5.2f | %.2e\n", loads(i), d_local(i), ...
            r_local(i), d_trace(i), r_trace(i), d_local(i)/scale(i));
    end
    n_failed = n_failed + check(all(r_local(2:end) >= 1.9) && all(r_trace(2:end) >= 1.9), ...
        material + ": difference of order load^2");
    n_failed = n_failed + check(all(d_local./scale <= 5*loads), ...
        material + ": relative difference below 5 load");
end
fprintf("\n%d check(s) failed.\n", n_failed);

function X = local_dofs(solution)
% [σ; u; γ] dofs, one column per element
dofs = arrayfun(@(f) as_matrix(f.dof), solution.Unknowns.field, 'UniformOutput', false);
X = cat(1, dofs{:});
end

function A = as_matrix(A)
if iscell(A), A = cat(1, A{:}); end
end
