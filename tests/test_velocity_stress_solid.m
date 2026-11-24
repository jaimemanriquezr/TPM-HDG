% test_velocity_stress_solid.m
% velocity_stress_solid_march (tests/helpers): the velocity-stress solid on
% the snapped carve of a moving disc (circle of radius 0.2, crossed grid),
% (sigma, v) transported by the exact two-mesh projection and advanced by
% one PEERS HDG trapezoid (or BDF2) solve on the new carve.
%   1. constitutive_source is option-gated: a solve with the field absent,
%      empty or zero gives the same dofs (bit for bit).
%   2. patch test: v and sigma polynomials of degree k in the body frame,
%      sigma linear in t, reproduced to round-off (1e-10 relative) at
%      k = 1, 2 on the fixed carve (N = 8, 16), on the translating carve
%      (A = 0.1, N = 16, 32) and on the translating and rotating carve
%      (Theta = 0.5, co-rotational transport), trapezoid and BDF2.
%   3. stability, homogeneous data and a random initial (sigma, v):
%      k = 1, N = 16, A = 0.1, dt = 0.025, 40 steps: no solve increases the
%      energy 1/2 (rho |v|^2 + (A sigma, sigma)) (1e-12 relative), no
%      transport increases it by more than 1e-3 relative, and the energy
%      decays over the run.
failures = 0;

% 1. option gate
bg = StructuredSimplicialGrid(8, 8, "crossed");
curve = CircularCurve([0.5; 0.5], 0.2);
[~, mesh] = bg.split(curve, max_normal_deviation = 90);
mesh = snap_to_curve(mesh, curve, normal_sign = -1);
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
mesh.faces(3, is_boundary) = mesh.boundary_tags("Dirichlet");
sp = LinearElasticityProblem(1, density = 1, time_step = 0.1, time_scheme = "newmark", ...
    young_modulus = 10, poisson_ratio = 0.3, source_term = {@(x, y) x.*y; @(x, y) x - y});
sp.Data.boundary_conditions = boundary_condition("Dirichlet", x_data = @(x, y) x, y_data = @(x, y) y.^2);
sp.Data.state = struct("solid", struct("displacement", ones(6, size(mesh.elements, 2))));
dofs = @(s) [reshape(s.Unknowns.field(1).dof, [], 1); reshape(s.Unknowns.field(2).dof, [], 1); ...
    reshape(s.Unknowns.hybrid(1).dof, [], 1)];
x0 = dofs(solve(sp, mesh, 1));
sp.Data.constitutive_source = [];
x1 = dofs(solve(sp, mesh, 1));
sp.Data.constitutive_source = zeros(14, size(mesh.elements, 2));
x2 = dofs(solve(sp, mesh, 1));
fprintf("1. gate: |absent - empty| %.1e, |absent - zero| %.1e\n", norm(x0 - x1), norm(x0 - x2));
failures = failures + check(isequal(x0, x1) && isequal(x0, x2), "constitutive_source gate");

% 2. patch test
cases = {struct("A", 0, "Theta", 0, "N", [8 16]), struct("A", 0.1, "Theta", 0, "N", [16 32]), ...
    struct("A", 0.1, "Theta", 0.5, "N", [16 32])};
for c = 1:numel(cases)
    for k = 1:2
        for N = cases{c}.N
            for scheme = ["trapezoid", "bdf2"]
                o = velocity_stress_solid_march(k, N, cases{c}.A, 0.05, 4, fields = "patch", ...
                    rotation = cases{c}.Theta, scheme = scheme);
                r = max(o.e_sigma/o.n_sigma, o.e_v/o.n_v);
                fprintf("2. patch A = %g Theta = %g k = %d N = %d %s: relative error %.2e\n", ...
                    cases{c}.A, cases{c}.Theta, k, N, scheme, r);
                failures = failures + check(r < 1e-10, sprintf("patch A = %g Theta = %g k = %d N = %d %s", ...
                    cases{c}.A, cases{c}.Theta, k, N, scheme));
            end
        end
    end
end

% 3. short stability check
o = velocity_stress_solid_march(1, 16, 0.1, 0.025, 40, homogeneous = true);
E = o.energy;
solve_rise = max((E(2:end) - o.energy_transported)./o.energy_transported);
transport_rise = max((o.energy_transported - E(1:end - 1))./E(1:end - 1));
fprintf("3. stability: max solve rise %.2e, max transport rise %.2e, E(end)/E(1) %.3e\n", ...
    solve_rise, transport_rise, E(end)/E(1));
failures = failures + check(solve_rise < 1e-12, "solve energy non-increasing");
failures = failures + check(transport_rise < 1e-3, "transport energy rise below 1e-3");
failures = failures + check(E(end) < E(1), "energy decays");

if failures == 0
    fprintf("test_velocity_stress_solid: all checks passed\n");
else
    fprintf("test_velocity_stress_solid: %d check(s) failed\n", failures);
end
