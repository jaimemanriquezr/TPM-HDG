% test_eulerian_hyperelastic_fsi.m
% eulerian_hyperelastic_fsi_march (tests/helpers): the fully Eulerian
% hyperelastic solid coupled to the PEERS Stokes fluid by
% FluidStructureCoupling.interface_problem (Data.eulerian_blocks), the
% circle of radius 0.2 translating with A = 0.1, both sides carved from
% the crossed grid, BDF2, Newton around the coupled solve; fluid data of
% the coupled velocity-stress test, solid fields "standard" with a bounded
% w; the displacement skeleton closed on Sigma with the kinematic g_u.
%   1. the Newton operator of the kinematic g_u: central differences of the
%      solid residual (eulerian_hyperelastic_solid_march's jacobian_check,
%      characteristic closure with data.kinematic, rotation 0.5, N = 8)
%      against the operator: k = 1 linear law (the residual is quadratic:
%      exact up to round-off) at most 1e-8 at every step (measured at most
%      3.0e-9); k = 2 neo-Hookean: at most 1e-7 at the step 1e-5 and
%      falling by at least 50 from 1e-3 to 1e-4 (measured 1.9e-8, 100).
%   2. linear law, k = 1, dt = 0.025, 3 steps, N = 16, 32: Newton converges
%      in at most 5 coupled solves per step (measured 4), and the errors at
%      T of u_F and v_S fall by at least 1.5 from N = 16 to 32 (measured
%      4.0 and 2.1).
%   3. k = 2, N = 16, dt = 0.0125, 8 steps: Newton in at most 5 solves per
%      step (measured 3-4); v_S error at T at most 3e-3 (measured 1.6e-3
%      with the full extension, 1.5e-3 with extension_degree = 1;
%      the extrapolated g_u gives 8.8e-3, the exact g_u 8.4e-4).
%   4. the neo-Hookean and Saint Venant-Kirchhoff laws (k = 1, N = 16,
%      3 steps): Newton converges.
%   5. flat interface (y = 0.45 + 0.1 sin(2 pi t), solid below, walls with
%      the exact u_hat, v_hat), k = 1, dt = 0.0125, 4 steps, N = 16, 32:
%      Newton in at most 5 solves per step (measured 4), the errors at T of
%      u_F and v_S fall by at least 1.5 (measured 3.2 and 1.7).
%   6. time order 2: data that start from rest (fluid and solid fields
%      times t^3), fixed carve (A = 0), k = 2, N = 16, T = 0.4, dt = 0.025,
%      0.0125, 0.00625: the self-convergence orders (dof 2-norm of the
%      difference of the final states) of u_F, v_S, u_S, sigma_S are at
%      least 1.8 (measured 2.01, 1.94, 1.93, 2.13). With data that do not
%      start from rest the projected initial levels excite a component of
%      the spatial error that depends on dt and does not decay over T = 1,
%      and the self-differences of v_S, sigma_S, u_F stall at that level.
%   7. perturbed twin (neo-Hookean, A = 0.1, k = 1, N = 16, dt = 0.025, 40
%      steps, initial_perturbation = [7, 1e-6], Newton tolerance 1e-12):
%      the solid's perturbation energy falls below 1e-2 of its start and
%      does not grow over the second half (it settles near 3e-4 of its
%      start, and stays there over 8 periods in longer runs).
%   8. Navier-Stokes fluid (convection = true: the fluid source gains rho_f
%      (u.grad) u, the fluid convected by the divergence-free postprocessed
%      velocity u* of the last iterate), SVK, k = 2, N = 8, 16, 32, dt =
%      0.2/N, T = 0.1: Newton (with the Picard update of u*) in at most 8
%      solves per step (measured 4-7), max_K ||div u*|| at most 1e-10
%      (measured 5e-13), and the EOC from N = 16 to 32 at least 1.9 for
%      u_F, v_S, u_S, sigma_S (measured 2.14, 2.02, 2.12, 2.03; the Stokes
%      run of the same data 2.43, 2.02, 2.12, 2.03: dt ~ h, so BDF2 caps
%      the combined order at 2).
failures = 0;
uu ="sin(pi*x)*sin(pi*y)*cos(t)";
pp = "cos(pi*x)*(y - 0.5)*(1 + t)";

% 1. Newton operator of the kinematic g_u
cases = {{1, "linear"}, {2, "neo_hookean"}};
for c = cases
    [k, law] = c{1}{:};
    dt = 0.025;
    ops = eulerian_hyperelastic_solid_march(k, 8, 0.1, dt, 0, rotation = 0.5, law = law, ...
        displacement_boundary = "characteristic", operators = true);
    ex = ops.exact;
    d2 = nchoosek(k + 2, 2);
    curve = CircularCurve([0.5; 0.5], 0.2);
    [~, sm] = StructuredSimplicialGrid(8, 8, "crossed").split(curve, max_normal_deviation = 90);
    sm = snap_to_curve(sm, curve, normal_sign = -1);
    geo = ops.geometry(sm);
    level0 = ops.project_vector(sm, k, @(x, y) [ex.v(x, y, dt); ex.u(x, y, dt)]);
    Z = ops.project_vector(sm, k, @(x, y) [ex.v(x, y, 0); ex.u(x, y, 0)]);
    data = struct("f", ops.moments(sm, k, @(x, y) ex.f(x, y, dt)), ...
        "g", ops.moments(sm, k, @(x, y) ex.g(x, y, dt)), "lam", zeros(4*geo.n_dof, 1), ...
        "dir", zeros(0, 1), "gu_face", zeros(size(geo.Pf, 1), size(sm.elements, 2), 2), ...
        "kinematic", struct("Zu", Z(2*d2 + 1:end, :), "dt", dt));
    lam0 = ops.initial_skeleton(geo, k, level0, data);
    je = ops.jacobian_check(geo, k, ops.hp, ops.prm, data, Z, level0, lam0);
    fprintf("1. kinematic g_u, k = %d %s: relative error at steps 1e-3..1e-8: %s\n", k, law, ...
        strjoin(compose("%.1e", je), " "));
    if law == "linear"
        failures = failures + check(max(je) <= 1e-8, "Jacobian, linear law");
    else
        failures = failures + check(je(3) <= 1e-7 && je(1)/je(2) >= 50, "Jacobian, neo-Hookean");
    end
end

% 2. k = 1, N = 16, 32
e = zeros(2, 4);
Ns = [16 32];
for i = 1:2
    o = eulerian_hyperelastic_fsi_march(1, Ns(i), 0.1, 0.025, 3, uu, pp);
    e(i, :) = o.step_errors(end, :);
    fprintf("2. N = %d: errors at T [u_F v_S u_S sigma_S] %s, Newton iterations %s\n", Ns(i), ...
        mat2str(e(i, :), 3), mat2str(o.iterations.'));
    failures = failures + check(all(o.converged) && max(o.iterations) <= 5, ...
        sprintf("Newton, N = %d", Ns(i)));
end
ratio = e(1, :)./e(2, :);
fprintf("2. error ratios N = 16 / 32: %s\n", mat2str(ratio, 3));
failures = failures + check(ratio(1) >= 1.5 && ratio(2) >= 1.5, "errors fall with h");

% 3. k = 2
o = eulerian_hyperelastic_fsi_march(2, 16, 0.1, 0.0125, 8, uu, pp);
fprintf("3. k = 2: errors at T %s, Newton iterations %s\n", mat2str(o.step_errors(end, :), 3), ...
    mat2str(o.iterations.'));
failures = failures + check(all(o.converged) && max(o.iterations) <= 5 && o.step_errors(end, 2) <= 3e-3, ...
    "k = 2 with the kinematic g_u");

% 4. neo-Hookean and Saint Venant-Kirchhoff
for law = ["neo_hookean", "svk"]
    o = eulerian_hyperelastic_fsi_march(1, 16, 0.1, 0.025, 3, uu, pp, law = law);
    fprintf("4. %s: errors at T %s, Newton iterations %s\n", law, mat2str(o.step_errors(end, :), 3), ...
        mat2str(o.iterations.'));
    failures = failures + check(all(o.converged) && all(isfinite(o.step_errors(:))), law);
end

% 5. flat interface
ef = zeros(2, 4);
for i = 1:2
    o = eulerian_hyperelastic_fsi_march(1, Ns(i), 0.1, 0.0125, 4, uu, pp, geometry = "flat");
    ef(i, :) = o.step_errors(end, :);
    fprintf("5. flat N = %d: errors at T %s, Newton iterations %s\n", Ns(i), mat2str(ef(i, :), 3), ...
        mat2str(o.iterations.'));
    failures = failures + check(all(o.converged) && max(o.iterations) <= 5, sprintf("flat Newton, N = %d", Ns(i)));
end
ratio = ef(1, :)./ef(2, :);
fprintf("5. flat error ratios N = 16 / 32: %s\n", mat2str(ratio, 3));
failures = failures + check(ratio(1) >= 1.5 && ratio(2) >= 1.5, "flat errors fall with h");

% 6. time order from rest
uu0 = "sin(pi*x)*sin(pi*y)*t^3";
pp0 = "cos(pi*x)*(y - 0.5)*t^3";
w0 = ["0.02*sin(5*X + 3*Y)*t^3"; "0.02*cos(4*X - 5*Y)*t^3"];
dts = [0.025 0.0125 0.00625];
r = cell(1, 3);
for i = 1:3
    r{i} = eulerian_hyperelastic_fsi_march(2, 16, 0, dts(i), round(0.4/dts(i)), uu0, pp0, w = w0);
end
d = zeros(2, 4);
for i = 1:2
    d(i, :) = [norm(r{i}.fluid_velocity - r{i + 1}.fluid_velocity, "fro"), norm(r{i}.v - r{i + 1}.v, "fro"), ...
        norm(r{i}.u - r{i + 1}.u, "fro"), norm(r{i}.sigma - r{i + 1}.sigma, "fro")];
end
order = log2(d(1, :)./d(2, :));
fprintf("6. from rest, fixed carve, k = 2, N = 16: self-convergence orders [u_F v_S u_S sigma_S] %s\n", ...
    mat2str(order, 3));
failures = failures + check(all(order >= 1.8), "time order 2 from rest");

% 7. perturbed twin
args = {"law", "neo_hookean", "tolerance", 1e-12};
r = eulerian_hyperelastic_fsi_march(1, 16, 0.1, 0.025, 40, uu, pp, args{:}, record_state = true);
o = eulerian_hyperelastic_fsi_march(1, 16, 0.1, 0.025, 40, uu, pp, args{:}, ...
    initial_perturbation = [7, 1e-6], reference = r.record);
P = [o.perturbation_energy; o.difference(:, 1)];
fprintf("7. twin: P(t)/P0 at t = 0.25, 0.5, 0.75, 1: %s, max fluid difference energy %.2e\n", ...
    mat2str(P([11 21 31 41]).'/P(1), 3), max(o.difference(:, 2)));
failures = failures + check(P(end) <= 1e-2*P(1) && max(P(21:end)) <= 1.05*P(21), "twin decays, then no growth");
% 8. Navier-Stokes fluid
e = zeros(3, 4);
Ns = [8 16 32];
for i = 1:3
    o = eulerian_hyperelastic_fsi_march(2, Ns(i), 0.1, 0.2/Ns(i), round(Ns(i)/2), uu, pp, law = "svk", ...
        convection = true);
    e(i, :) = o.step_errors(end, :);
    fprintf("8. Navier-Stokes N = %d: errors at T %s, Newton %s, max div u* %.1e\n", Ns(i), ...
        mat2str(e(i, :), 3), mat2str(o.iterations.'), max(o.divergence_star));
    failures = failures + check(all(o.converged) && max(o.iterations) <= 8 && max(o.divergence_star) <= 1e-10, ...
        sprintf("Navier-Stokes Newton, N = %d", Ns(i)));
end
order = log2(e(2, :)./e(3, :));
fprintf("8. Navier-Stokes EOC (N = 16 / 32) %s\n", mat2str(order, 3));
failures = failures + check(all(order >= 1.9), "Navier-Stokes coupled order 2");

if failures == 0
    fprintf("test_eulerian_hyperelastic_fsi: all checks passed\n");
else
    fprintf("test_eulerian_hyperelastic_fsi: %d check(s) failed\n", failures);
end
