function pde_elasticity_mixed_bc(k, max_iterations)
% Component-wise boundary-condition verification for linear elasticity.
%
% On the right side x = 1 (outward normal n = (1,0)) the trace is fixed
% component-wise: u_x = exact (dirichlet) and (sigma n)_y = sigma_21 exact
% (neumann); all other boundaries carry full Dirichlet data. The errors
% must match the all-Dirichlet baseline and converge at the same rate.
% Sign convention verified here: neumann data is the outward-normal
% traction component (sigma n)_c, positive sign.
arguments
    k               = 2
    max_iterations  = 4
end
u_exact = ["exp(x)*sin(y)"; "exp(y)*cos(x)"];   % sigma_21 nonzero on x = 1
E = 1e2; nu = 0.3;
data = LinearElasticityProblem.get_manufactured_solution(u_exact, E, nu);
problem = LinearElasticityProblem(young_modulus = E, poisson_ratio = nu, ...
    source_term = data.f, datum_dirichlet = data.u);

[e_dir, e_mix, h_list] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    N = 2^(i + 1);
    mesh = StructuredSimplicialGrid(N, N);
    mesh = apply_tag(mesh, "Mixed", @(X) isapprox(X(1, :), 1));
    h_list(i) = sqrt(2)/N;

    pb = problem;
    pb.Data.boundary_conditions = [ ...
        boundary_condition("Dirichlet", x_data = data.u{1}, y_data = data.u{2});
        boundary_condition("Mixed",     x_data = data.u{1}, y_data = data.u{2})];
    sol = solve(pb, mesh, Order= k);
    err = sol.compute_l2_error({[], data.u});
    e_dir(i) = err(2);

    pm = problem;
    pm.Data.boundary_conditions = [ ...
        boundary_condition("Dirichlet", x_data = data.u{1}, y_data = data.u{2});
        boundary_condition("Mixed", y = "neumann", x_data = data.u{1}, ...
            y_data = data.sigma{2, 1})];
    sol = solve(pm, mesh, Order= k);
    err = sol.compute_l2_error({[], data.u});
    e_mix(i) = err(2);
end

r_mix = [NaN; compute_eoc(e_mix, h_list)];
fprintf("\n%6s %12s %12s %8s %8s\n", "N", "e_u(dir)", "e_u(mixed)", "eoc", "ratio");
for i = 1:max_iterations
    fprintf("%6d %12.4e %12.4e %8.2f %8.3f\n", ...
        2^(i + 1), e_dir(i), e_mix(i), r_mix(i), e_mix(i)/e_dir(i));
end
fprintf("PASS requires: ratio ~ 1 on every row and eoc ~ %d.\n", k + 1);
end
