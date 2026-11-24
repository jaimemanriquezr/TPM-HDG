% test_eulerian_hyperelastic_solid.m
% eulerian_hyperelastic_solid_march (tests/helpers): the fully Eulerian
% hyperelastic solid in (v, u) on the snapped carve of a disc (circle of
% radius 0.2, crossed grid), advective HDG transport, BDF2, the old levels
% carried by the zero-shift two-mesh projection.
%   1. patch test: u, v polynomials in x (u (1 + t^2), v (1 + t)), so that
%      the exact solution is in the discrete space; relative L2 errors of
%      v, u and sigma at the final time (each divided by the L2 norm of the
%      exact field on the final carve) at most 1e-10:
%      a. fixed carve, linear law, degree-k fields, N = 8, 16,
%         displacement_boundary "dirichlet" (k = 1, 2, 3) and "free"
%         (k = 1, 2; at k = 3 the free boundary, u_hat the inner trace on
%         the inflow faces of the transport of u, amplifies round-off by a
%         factor of about 30 per step at N = 16, dt = 0.05);
%      b. fixed carve, neo-Hookean law, degree-1 fields, k = 1, 2, 3, N = 8, 16;
%      c. moving carve (A = 0.1; A = 0.1 with rotation 0.5), degree-1
%         fields with extension_degree = 1, linear and neo-Hookean,
%         k = 1, 2, 3, N = 16, 32;
%      d. moving carve (A = 0.1), degree-k fields with the full extension
%         (extension_degree = []), linear law, k = 2, 3, N = 16, 32.
%      dt = 0.05, 3 steps, picard_tolerance = 1e-12 (Newton; the lagged
%      Picard loop converges only linearly, at a rate of up to about 0.4
%      on the neo-Hookean cases at N = 32: stopped at a relative change of
%      1e-10, its own error reaches 3e-10 in sigma there).
%   2. stability: k = 1, N = 16, A = 0.1, dt = 0.025, 40 steps, neo-Hookean,
%      rigid motion (w = 0), a perturbed twin of linearized energy 1e-12
%      (initial_perturbation = [7, 1e-6]): every Picard loop converges and
%      the perturbation energy at the end is below its initial value.
%   3. characteristic closure of u_hat on the boundary: patch fields over
%      12 steps at N = 32 stay at round-off (k = 1, 3).
%   4. Newton on the transport: (a) jacobian_check, "quadratic" fields,
%      neo-Hookean, N = 8, k = 1, 2: the central difference of the
%      residual against the Newton operator has relative error at most
%      1e-6 at the step 1e-4 and falls by at least 50 from the step 1e-3 to
%      1e-4 (second order: the derivative is exact; measured 1.4e-7 and
%      6.0e-7, ratio 100); (b) rigid translation and rotation (A = 0.1,
%      rotation 0.5, w = 0), k = 1, N = 16, dt = 0.025, 12 steps: Newton
%      converges in at most 5 iterations per step (measured 4; the lagged
%      Picard loop needs 10 to 30, its count growing with |grad u|).
%   5. initial iterate under rotation: k = 3, neo-Hookean, A = 0.1,
%      rotation 0.5, w = 0, N = 16, dt = 0.025, 4 steps run and converge
%      (with the data of t_n+1 as the boundary u_hat of the initial
%      iterate, det(I - L_h) reached -0.1 there and the first assembly
%      failed).
%   6. the Saint Venant-Kirchhoff law: patch test and Jacobian (below).
failures = 0;
dt = 0.05;
n_steps = 3;
bound = 1e-10;
report = @(o) [o.e_v/o.n_v, o.e_u/o.n_u, o.e_sigma/o.n_sigma];

% 1. patch test
cases = struct("name", {}, "A", {}, "Theta", {}, "law", {}, "pdeg", {}, "ext", {}, ...
    "boundary", {}, "k", {}, "N", {});
cases(end + 1) = struct("name", "a. fixed, linear", "A", 0, "Theta", 0, "law", "linear", ...
    "pdeg", [], "ext", 1, "boundary", "dirichlet", "k", 1:3, "N", [8 16]);
cases(end + 1) = struct("name", "a. fixed, linear", "A", 0, "Theta", 0, "law", "linear", ...
    "pdeg", [], "ext", 1, "boundary", "free", "k", 1:2, "N", [8 16]);
cases(end + 1) = struct("name", "b. fixed, neo-Hookean", "A", 0, "Theta", 0, "law", "neo_hookean", ...
    "pdeg", 1, "ext", 1, "boundary", "dirichlet", "k", 1:3, "N", [8 16]);
for Theta = [0, 0.5]
    for law = ["linear", "neo_hookean"]
        cases(end + 1) = struct("name", "c. moving, degree 1", "A", 0.1, "Theta", Theta, "law", law, ...
            "pdeg", 1, "ext", 1, "boundary", "dirichlet", "k", 1:3, "N", [16 32]); %#ok<SAGROW>
    end
end
cases(end + 1) = struct("name", "d. moving, degree k, full extension", "A", 0.1, "Theta", 0, ...
    "law", "linear", "pdeg", [], "ext", [], "boundary", "dirichlet", "k", 2:3, "N", [16 32]);
for c = cases
    for k = c.k
        for N = c.N
            o = eulerian_hyperelastic_solid_march(k, N, c.A, dt, n_steps, rotation = c.Theta, ...
                law = c.law, fields = "patch", patch_degree = c.pdeg, extension_degree = c.ext, ...
                displacement_boundary = c.boundary, picard_tolerance = 1e-12);
            r = report(o);
            label = sprintf("%s A = %g Theta = %g %s %s k = %d N = %d", c.name, c.A, c.Theta, ...
                c.law, c.boundary, k, N);
            fprintf("1%s: relative errors v %.2e u %.2e sigma %.2e, Picard %s\n", label, r, ...
                mat2str(o.iterations.'));
            failures = failures + check(max(r) <= bound && all(o.converged), label);
        end
    end
end

% 2. stability of a perturbation of the rigid motion
o = eulerian_hyperelastic_solid_march(1, 16, 0.1, 0.025, 40, law = "neo_hookean", ...
    w = ["0"; "0"], initial_perturbation = [7, 1e-6]);
P = o.perturbation_energy;
fprintf("2. stability: perturbation energy initial %.3e final %.3e (ratio %.3e), max ratio %.3e, " + ...
    "Picard iterations %d..%d, all converged %d\n", P(1), P(end), P(end)/P(1), max(P)/P(1), ...
    min(o.iterations), max(o.iterations), all(o.converged));
failures = failures + check(all(o.converged), "Picard converges");
failures = failures + check(P(end) < P(1), "perturbation energy decays");

% 3. characteristic closure of u_hat on the boundary (no displacement data
% on outflow): patch fields, linear law, fixed carve, N = 32, dt = 0.05,
% 12 steps (the "free" closure grows to a relative sigma error of order 1
% here at k = 1 and 5e4 at k = 3); relative errors at most 1e-10 at the end
for k = [1 3]
    o = eulerian_hyperelastic_solid_march(k, 32, 0, 0.05, 12, law = "linear", fields = "patch", ...
        displacement_boundary = "characteristic", picard_tolerance = 1e-12);
    r = report(o);
    label = sprintf("3. characteristic k = %d N = 32, 12 steps", k);
    fprintf("%s: relative errors v %.2e u %.2e sigma %.2e, iterations %d..%d\n", label, r, ...
        min(o.iterations), max(o.iterations));
    failures = failures + check(max(r) <= bound && all(o.converged), label);
end

% 4. Newton on the transport
for k = 1:2
    o = eulerian_hyperelastic_solid_march(k, 8, 0, 0.05, 1, fields = "quadratic", ...
        jacobian_check = true, picard_max_iterations = 1);
    je = o.jacobian_error;
    fprintf("4a. Jacobian k = %d: relative error at steps 1e-3..1e-8: %s\n", k, ...
        strjoin(compose("%.1e", je), " "));
    failures = failures + check(je(2) <= 1e-6 && je(1)/je(2) >= 50, sprintf("Jacobian k = %d", k));
end
o = eulerian_hyperelastic_solid_march(1, 16, 0.1, 0.025, 12, rotation = 0.5, w = ["0"; "0"]);
fprintf("4b. Newton under rotation: iterations %s, all converged %d\n", mat2str(o.iterations.'), ...
    all(o.converged));
failures = failures + check(all(o.converged) && max(o.iterations) <= 5, "Newton iterations under rotation");

% 5. initial iterate under rotation, k = 3
try
    o = eulerian_hyperelastic_solid_march(3, 16, 0.1, 0.025, 4, rotation = 0.5, w = ["0"; "0"]);
    fprintf("5. k = 3 rotation: iterations %s, all converged %d, relative errors v %.2e u %.2e\n", ...
        mat2str(o.iterations.'), all(o.converged), o.e_v/o.n_v, o.e_u/o.n_u);
    failures = failures + check(all(o.converged), "k = 3 under rotation");
catch err
    fprintf("5. k = 3 rotation: %s\n", err.message);
    failures = failures + check(false, "k = 3 under rotation");
end

% 6. Saint Venant-Kirchhoff (law = "svk"): (a) patch test on the moving
% carve (A = 0.1, rotation 0.5), degree-1 fields, extension_degree = 1,
% k = 1, 2, 3, N = 8, 16, relative errors at most 1e-10 (measured at most
% 3e-13); (b) Jacobian check as in 4a, at most 1e-7 at the step 1e-5
% (measured 1.8e-8 at k = 1, 8.8e-8 at k = 2)
for k = 1:3
    for N = [8 16]
        o = eulerian_hyperelastic_solid_march(k, N, 0.1, dt, n_steps, rotation = 0.5, law = "svk", ...
            fields = "patch", patch_degree = 1, extension_degree = 1, picard_tolerance = 1e-12);
        r = report(o);
        label = sprintf("6a. SVK patch, moving, k = %d N = %d", k, N);
        fprintf("%s: relative errors v %.2e u %.2e sigma %.2e\n", label, r);
        failures = failures + check(max(r) <= bound && all(o.converged), label);
    end
end
for k = 1:2
    o = eulerian_hyperelastic_solid_march(k, 8, 0, 0.05, 1, fields = "quadratic", law = "svk", ...
        jacobian_check = true, picard_max_iterations = 1);
    je = o.jacobian_error;
    fprintf("6b. SVK Jacobian k = %d: relative error at steps 1e-3..1e-8: %s\n", k, ...
        strjoin(compose("%.1e", je), " "));
    failures = failures + check(je(3) <= 1e-7, sprintf("SVK Jacobian k = %d", k));
end

if failures == 0
    fprintf("test_eulerian_hyperelastic_solid: all checks passed\n");
else
    fprintf("test_eulerian_hyperelastic_solid: %d check(s) failed\n", failures);
end
