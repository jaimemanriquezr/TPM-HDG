% test_velocity_stress_fsi.m
% velocity_stress_fsi_march (tests/helpers): the FSI march with the
% velocity-stress Eulerian solid, fluid and solid carved from one
% background grid at every step, Sigma(t) moving rigidly, the solid state
% (sigma, v) transported exactly and coupled to the PEERS fluid by
% FluidStructureCoupling.interface_problem with the solid's
% constitutive_source; schemes "midpoint" (Crank-Nicolson at the midpoint
% on both sides) and "bdf2" (first step midpoint).
%   1. the constitutive_source gate of interface_problem: one coupled step
%      (flat, "left" N = 8, k = 1) with the field absent, empty and zero
%      gives the same dofs bit for bit (measured difference 0); a random
%      nonzero one changes them (by 893); the solid problem interface_problem
%      builds, solved alone on an all-Dirichlet carve, agrees with the
%      low-level LinearElasticityProblem solve with the same source to
%      4.1e-15 (relative, stress and velocity).
%   2. patch test, both schemes: fluid stream function of degree k + 1 and
%      pressure of degree k, solid fields = "patch" (body-frame V of degree
%      k, Sigma = S0 + t C eps(V)), all linear in t; dt = 0.05, 9 steps
%      (Sigma crosses y = 0.5 up and back: rows switch on and off on both
%      sides). Flat, "left" N = 8, 16, A = 0 and 0.1, k = 1, 2 (and N = 8,
%      k = 3); circle,
%      "crossed" N = 16, 32 (k = 1) and N = 16 (k = 2), A = 0, A = 0.1 and
%      A = 0.1 with Theta = 0.5. The worst step error of sigma_F, u_F,
%      sigma_S, v_S relative to the exact norm is below 1e-10 (measured at
%      most 1.5e-12 on these cases, the curved Sigma included: the
%      transfer along the paths is exact for polynomials, so the error is
%      round-off).
%   3. stability, circle, crossed N = 16, k = 1, A = 0.1, dt = 0.025, 40
%      steps, both schemes: (i) smooth data, the worst step error of each
%      field is at most 3 times its worst over the first 4 steps (measured
%      (sigma_F, u_F, sigma_S, v_S): midpoint [1.06 1 1.38 1.39], bdf2
%      [1 1 1.73 1]; the exact fields themselves change over T = 1);
%      (ii) the march with initial_perturbation = [1, 1e-6] against the
%      unperturbed one (record_state, reference): the solid energy of the
%      difference at the end is not larger than at the start (measured
%      end/start 3.65e-3 midpoint, 3.83e-3 bdf2; its largest value over
%      the steps 0.195 of the start for both).
%   4. (A = 0 against one fixed carve) skipped: at A = 0 the march executes
%      no transport and the fluid remap copies every element's dofs, so
%      there is no executed-with-identity path to compare.
%   5. restart: 4 + 6 steps (start_step, start_state) reproduce 10 steps
%      (circle, N = 16, k = 1, A = 0.1, Theta = 0.5, dt = 0.025) exactly:
%      isequal on the final state, the step errors, the energies and the
%      records, and a zero difference to the uninterrupted run's record.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
uu = "sin(pi*x)*sin(pi*y)*cos(t)";  pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
us = ["0.1*sin(2*x + y)*cos(t)"; "0.1*cos(x - y)*(1 + t^2)"];

%% 1. constitutive_source gate
k = 1;  dt = 0.05;
[fp, ~, src, ~, curve] = translating_flat_problems(k, 0.45, 0, dt, dt, "x*y*(1 + t)", ...
    "(x + y)*(1 + t)", ["x*y"; "x - y"], 1, 1, 1, 10, 0.3);
fp.Data.time_scheme = "euler";
[fm, sm] = StructuredSimplicialGrid(8, 8, "left").split(curve, max_normal_deviation = 90);
fm = snap_to_curve(fm, curve);
sm = snap_to_curve(sm, curve, normal_sign = -1);
n_el = size(sm.elements, 2);
nP = PEERSBasis(k).Dimension;
stream = RandStream("mt19937ar", Seed = 3);
sp = LinearElasticityProblem(k, density = 1.5, time_step = dt, time_scheme = "newmark", ...
    newmark_beta = 1, newmark_gamma = 0.5, young_modulus = 10/1.5, poisson_ratio = 0.3, ...
    source_term = {@(x, y) x.*y; @(x, y) x - y}, datum_dirichlet = {@(x, y) x; @(x, y) y.^2});
state = struct("fluid", struct("velocity", zeros(6, size(fm.elements, 2)), "acceleration", 0), ...
    "solid", struct("displacement", randn(stream, 6, n_el), "velocity", 0, "acceleration", 0, ...
    "stress", randn(stream, nP, n_el), "rotation", 0));
x0 = coupled_dofs(fp, sp, fm, sm, curve, src, state, k);
sp.Data.constitutive_source = [];
x1 = coupled_dofs(fp, sp, fm, sm, curve, src, state, k);
sp.Data.constitutive_source = zeros(nP, n_el);
x2 = coupled_dofs(fp, sp, fm, sm, curve, src, state, k);
sp.Data.constitutive_source = randn(stream, nP, n_el);
x3 = coupled_dofs(fp, sp, fm, sm, curve, src, state, k);
fprintf("1. gate: |absent - empty| %.1e, |absent - zero| %.1e, |nonzero - absent| %.3g\n", ...
    max(abs(x0 - x1)), max(abs(x0 - x2)), max(abs(x3 - x0)));
n_fail = n_fail + check(isequal(x0, x1) && isequal(x0, x2), "constitutive_source gate");
n_fail = n_fail + check(max(abs(x3 - x0)) > 1e-3, "a nonzero constitutive_source changes the solution");
is_boundary = accumarray(sm.boundaries(:), 1, [size(sm.faces, 2), 1]) == 1;
sd = sm;
sd.faces(3, is_boundary) = sd.boundary_tags("Dirichlet");
sp.Data.boundary_conditions = boundary_condition("Dirichlet", x_data = @(x, y) x, y_data = @(x, y) y.^2);
[~, ~, solid_pb] = FluidStructureCoupling.interface_problem(fp, sp, fm, sd, curve, src{:}, "onestep", ...
    "state", state);
declarative = solve(solid_pb, sd, k);
sp.Data.state = struct("solid", state.solid);
low = solve(sp, sd, k);
d_low = 0;
for i = 1:2
    a = declarative.Unknowns.field(i).dof;
    b = low.Unknowns.field(i).dof;
    if iscell(b), b = cat(1, b{:}); end
    d_low = max(d_low, max(abs(a(:) - b(:)))/max(abs(b(:))));
end
fprintf("1. solid alone, declarative vs low level: relative difference %.2e\n", d_low);
n_fail = n_fail + check(d_low < 1e-12, "declarative constitutive_source differs from the low-level one");

%% 2. patch test
cases = {"flat", 0, 0, [8 16], 1:2; "flat", 0.1, 0, [8 16], 1:2; ...
    "circle", 0, 0, [16 32], 1; "circle", 0.1, 0, [16 32], 1; "circle", 0.1, 0.5, [16 32], 1; ...
    "circle", 0, 0, 16, 2; "circle", 0.1, 0, 16, 2; "circle", 0.1, 0.5, 16, 2; ...
    "flat", 0, 0, 8, 3; "flat", 0.1, 0, 8, 3};
for scheme = ["midpoint", "bdf2"]
    for c = 1:size(cases, 1)
        [geometry, A, Theta, Ns, ks] = cases{c, :};
        style = "crossed";
        if geometry == "flat", style = "left"; end
        for k = ks
            [uu_p, pp_p] = patch_fluid(k);
            for N = Ns
                o = velocity_stress_fsi_march(k, StructuredSimplicialGrid(N, N, style), A, 0.05, 9, ...
                    uu_p, pp_p, us, scheme = scheme, geometry = geometry, rotation = Theta, ...
                    fields = "patch");
                rel = max(o.step_errors(:, 1:4), [], 1)./[o.n_sigma_F, o.n_u_F, o.n_sigma_S, o.n_v_S];
                name = sprintf("%s %s A = %g Theta = %g k = %d N = %d", scheme, geometry, A, Theta, k, N);
                fprintf("2. patch %s: relative %s, on/off %s\n", name, mat2str(rel, 2), ...
                    mat2str(sum(o.history(:, 4:5), 1)));
                n_fail = n_fail + check(max(rel) < 1e-10, "patch " + name);
                if A ~= 0
                    n_fail = n_fail + check(all(sum(o.history(:, 4:5), 1) > 0), ...
                        "patch " + name + ": no fluid element switched on and off");
                end
            end
        end
    end
end

%% 3. stability
bg = StructuredSimplicialGrid(16, 16, "crossed");
[factor_i, factor_ii] = deal(3, 1);
for scheme = ["midpoint", "bdf2"]
    run = @(varargin) velocity_stress_fsi_march(1, bg, 0.1, 0.025, 40, uu, pp, us, "scheme", scheme, ...
        "geometry", "circle", varargin{:});
    o = run("record_state", true);
    p = run("initial_perturbation", [1, 1e-6], "reference", o.record);
    growth = max(o.step_errors(:, 1:4), [], 1)./max(o.step_errors(1:4, 1:4), [], 1);
    % the perturbation's solid energy is 1e-12 times the initial energy
    Eg = p.difference(end, 1)/(1e-12*o.energy(1));
    fprintf("3. %s: worst / first-4-steps worst error (sF, uF, sS, vS) %s; difference solid energy " + ...
        "end/start %.3e (max over steps %.3e), fluid kinetic energy of the difference at the end %.3e\n", ...
        scheme, mat2str(growth, 3), Eg, max(p.difference(:, 1))/(1e-12*o.energy(1)), p.difference(end, 2));
    n_fail = n_fail + check(all(growth <= factor_i), scheme + ": smooth errors grow");
    n_fail = n_fail + check(Eg <= factor_ii, scheme + ": the perturbation grows");
end

%% 5. restart
for scheme = ["midpoint", "bdf2"]
    run = @(n, varargin) velocity_stress_fsi_march(1, bg, 0.1, 0.025, n, uu, pp, us, "scheme", scheme, ...
        "geometry", "circle", "rotation", 0.5, varargin{:});
    full = run(10, "record_state", true);
    a = run(4);
    b = run(10, "start_step", 4, "start_state", a.state, "record_state", true, "reference", full.record);
    same = isequal(full.state, b.state) && isequal(full.step_errors(5:end, :), b.step_errors) ...
        && isequal(full.energy(5:end), b.energy) && isequal(full.record(5:end), b.record) ...
        && all(b.difference == 0, "all");
    fprintf("5. restart %s: 4 + 6 steps isequal 10 steps: %d\n", scheme, same);
    n_fail = n_fail + check(same, scheme + ": the restart differs from the uninterrupted run");
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_velocity_stress_fsi:failed", "%d check(s) failed", n_fail);
end

function x = coupled_dofs(fp, sp, fm, sm, curve, src, state, k)
% all dofs of one coupled step
coupled = FluidStructureCoupling.interface_problem(fp, sp, fm, sm, curve, src{:}, "onestep", ...
    "Order", k, "state", state, "patch_weight", "jacobian", "swap_transmission", true);
[fs, ss] = solve(coupled, k);
c = [{fs.Unknowns.field.dof}, {fs.Unknowns.hybrid.dof}, {ss.Unknowns.field.dof}, {ss.Unknowns.hybrid.dof}];
c = cellfun(@(a) reshape(double(a), [], 1), c, UniformOutput = false);
x = cat(1, c{:});
end

function [uu, pp] = patch_fluid(k)
% stream function of degree k + 1 and pressure of degree k, linear in t
switch k
    case 1, uu = "(x^2 + x*y - 2*y^2 + x)*(1 + t)"; pp = "(x - 2*y + 1)*(1 + t)";
    case 2, uu = "(x^3 - 2*x*y^2 + x^2*y + y)*(1 + t)"; pp = "(x^2 - x*y + y)*(1 + t)";
    case 3, uu = "(x^4 - 3*x^2*y^2 + x*y^3 + y)*(1 + t)"; pp = "(x^3 - x*y^2 + y)*(1 + t)";
end
end
