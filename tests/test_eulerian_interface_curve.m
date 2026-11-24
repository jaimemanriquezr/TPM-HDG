% test_eulerian_interface_curve.m
% eulerian_interface_curve (tests/helpers): the current interface of an
% Eulerian solid, Sigma = {x : x - u(x) in Sigma_0}, as a DeformedCurve
% of the reference circle Sigma_0 (radius 0.2 about (0.5, 0.5)).
%   1. convergence in h: the "standard" fields of
%      eulerian_hyperelastic_solid_march (A = 0.1, rotation 0.5, the
%      default w, which does not vanish on Sigma_0) at t = 0.3; the exact
%      image curve from the exact u (samples every 0.002); the solid
%      carved against it (crossed grid, N = 16, 32, 64, snapped), the
%      L2 projection of the exact u (k = 2), and the curve rebuilt from
%      it: the largest distance between the two curves at equal reference
%      parameter falls at order at least 2.5 (measured 2.90, 3.28;
%      errors 2.7e-5, 3.7e-6, 3.8e-7).
%   2. the coupled circle march (eulerian_hyperelastic_fsi_march, SVK,
%      k = 2, N = 16, dt = 0.0125, 8 steps) with Sigma taken from the
%      solid's u^n, lagged one step (interface = "solid"): Newton
%      converges in at most 5 solves per step and the errors at T stay
%      within 1.5 times those of the prescribed circle (measured 0.76-0.85
%      times: the exact Sigma of these fields is not the circle).
%   3. a full period (k = 1, N = 16, dt = 0.025, 40 steps), interface =
%      "solid": every step converges and the largest error over the run
%      stays within 1.5 times the prescribed run's (measured [0.98 0.88
%      1.08 1.11]).
failures = 0;

% 1. convergence in h
k = 2;
ops = eulerian_hyperelastic_solid_march(k, 16, 0.1, 0.01, 0, rotation = 0.5, law = "svk", operators = true);
ex = ops.exact;
t = 0.3;
exact_u = @(x) reshape(ex.u(x(1, :).', x(2, :).', t), [], 2).';
reference = CircularCurve([0.5; 0.5], 0.2);
exact_curve = eulerian_interface_curve(reference, [], [], exact_u, spacing = 0.002);
ts = linspace(reference.t_range(1), reference.t_range(2), 2001);
Ns = [16 32 64];
e = zeros(size(Ns));
for i = 1:numel(Ns)
    [~, sm] = StructuredSimplicialGrid(Ns(i), Ns(i), "crossed").split(exact_curve, max_normal_deviation = 90);
    sm = snap_to_curve(sm, exact_curve, normal_sign = -1);
    U = ops.project_vector(sm, k, @(x, y) ex.u(x, y, t));
    [curve, info] = eulerian_interface_curve(reference, sm, k, U);
    e(i) = max(vecnorm(curve.get_2d_points(ts) - exact_curve.get_2d_points(ts)));
    fprintf("1. N = %d: max |x_h - x| = %.3e (fixed-point residual %.1e)\n", Ns(i), e(i), info.residual);
end
order = log2(e(1:end - 1)./e(2:end));
fprintf("1. orders %s\n", mat2str(order, 3));
failures = failures + check(all(order >= 2.5), "reconstructed Sigma converges in h");

% 2. coupled, k = 2
uu = "sin(pi*x)*sin(pi*y)*cos(t)";
pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
o = cell(1, 2);
interfaces = ["prescribed", "solid"];
for i = 1:2
    o{i} = eulerian_hyperelastic_fsi_march(2, 16, 0.1, 0.0125, 8, uu, pp, law = "svk", interface = interfaces(i));
    fprintf("2. %s: errors at T [u_F v_S u_S sigma_S] %s, Newton %s\n", interfaces(i), ...
        mat2str(o{i}.step_errors(end, :), 3), mat2str(o{i}.iterations.'));
end
ratio = o{2}.step_errors(end, :)./o{1}.step_errors(end, :);
fprintf("2. ratios solid / prescribed %s\n", mat2str(ratio, 3));
failures = failures + check(all(o{2}.converged) && max(o{2}.iterations) <= 5 && all(ratio <= 1.5), ...
    "coupled march with Sigma from the solid, k = 2");

% 3. a full period, k = 1
for i = 1:2
    o{i} = eulerian_hyperelastic_fsi_march(1, 16, 0.1, 0.025, 40, uu, pp, law = "svk", interface = interfaces(i));
end
ratio = max(o{2}.step_errors)./max(o{1}.step_errors);
fprintf("3. a period: largest errors, solid / prescribed %s, all converged %d\n", mat2str(ratio, 3), ...
    all(o{2}.converged));
failures = failures + check(all(o{2}.converged) && all(ratio <= 1.5), "a period with Sigma from the solid");

if failures == 0
    fprintf("test_eulerian_interface_curve: all checks passed\n");
else
    fprintf("test_eulerian_interface_curve: %d check(s) failed\n", failures);
end
