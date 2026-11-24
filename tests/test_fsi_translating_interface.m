% Transient FSI with a prescribed translating interface on fixed meshes:
% Sigma(t) = Sigma_0 + d(t), d(t) = A sin(2 pi t) e_x, a circle inside the
% gap of a fitted_split pair (both meshes retreated by `gap`, A < gap).
% Each step rebuilds both partitions against Sigma(t_n) (the curve passed
% to FluidStructureCoupling.solve) and translates the solid Newmark
% history by d(t_n) - d(t_n-1) (LinearElasticityProblem.translate_state);
% the data are translating_circle_problems'.
%
%   1. A = 0 reproduces the fixed-Sigma march (load_problems, the curve
%      Sigma_0, no translation) dof for dof (1e-12).
%   2. translate_state is exact on polynomials of degree k: the volume,
%      skeleton, stress and rotation translates equal the projections of
%      the translated field (1e-10).
%      translated_field is the translate of the volume levels: on the
%      rotation it equals translate_state's (1e-14), on the displacement
%      it is exact (1e-10).
%   3. The coupled march is exact on a solution both integrators and the
%      TPM transfer reproduce, for A > 0: u_S = p(x - d(t)) g(t) with p of
%      degree k and g quadratic (Newmark is exact along the
%      characteristics), u_F linear in t (backward Euler, then
%      Crank-Nicolson). Error below 1e-9 at every step, so the interface
%      data, the per-step partitions and the translated history are
%      consistent.
%
%See also: test_fsi_transient_polynomial, translating_circle_problems,
%          LinearElasticityProblem.translate_state
clc; clear;

k = 2; h = 0.1; gap = 0.05;
rho_f = 1; mu = 1; rho_s = 1; E = 10; nu = 0.3;
curve_0 = CircularCurve([0.5; 0.5], 0.4);
[fluid_mesh, solid_mesh] = fitted_split(curve_0, h, gap);
n_fail = 0;

%% 1. A = 0 against the fixed-Sigma march
uu = "sin(pi*x)*sin(pi*y)*cos(t)";
pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
us = ["x*sin(y)*cos(t)"; "y*cos(x)*(1 + t^2)"];
dt = 0.05; n_steps = 3;
ref = march(@(t) fixed_data(k, curve_0, t, dt, uu, pp, us, rho_f, mu, rho_s, E, nu), ...
    fluid_mesh, solid_mesh, k, dt, n_steps, []);
moving = march(@(t) translating_circle_problems(k, curve_0, 0, t, dt, uu, pp, us, ...
    rho_f, mu, rho_s, E, nu), fluid_mesh, solid_mesh, k, dt, n_steps, @(t) 0);
d_fix = max(cellfun(@(a, b) max(abs(a - b), [], "all"), ref.dofs, moving.dofs));
fprintf("1. A = 0 vs fixed Sigma: max dof difference %.2e\n", d_fix);
n_fail = n_fail + check(d_fix <= 1e-12, "A = 0 reproduces the fixed-Sigma march");

%% 2. translate_state on polynomials of degree k
solid = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, ...
    time_scheme = "newmark", density = rho_s);
shift = [0.004; -0.0025];
p = {@(x, y) x.^2 - 2*x.*y + y; @(x, y) y.^2 + x.*y - x};
p_shifted = cellfun(@(f) @(x, y) f(x - shift(1), y - shift(2)), p, UniformOutput = false);
fields_of = @(f) project(solid, solid_mesh, k, f);
state.solid = fields_of(p);
expected = fields_of(p_shifted);
state = translate_state(solid, solid_mesh, state, shift, k);
names = ["displacement", "displacement_skeleton", "stress", "rotation"];
d_poly = max(arrayfun(@(n) max(abs(state.solid.(n) - expected.(n)), [], "all"), names));
fprintf("2. translate_state on P_%d: max difference %.2e\n", k, d_poly);
n_fail = n_fail + check(d_poly <= 1e-10, "translate_state is exact on polynomials");
untranslated = fields_of(p);
gamma = @(g) cat(1, g{:});
if ~iscell(untranslated.rotation), gamma = @(g) g; end
d_rot = max(abs(translated_field(solid, solid_mesh, k, gamma(untranslated.rotation), shift) ...
    - gamma(state.solid.rotation)), [], "all");
d_disp = max(abs(translated_field(solid, solid_mesh, k, untranslated.displacement, shift) ...
    - expected.displacement), [], "all");
fprintf("2. translated_field: vs translate_state's rotation %.2e, on P_%d %.2e\n", d_rot, k, d_disp);
n_fail = n_fail + check(d_rot <= 1e-14 && d_disp <= 1e-10, ...
    "translated_field is translate_state's volume translate");

%% 3. coupled march, exact solution, A > 0
A = 0.01; dt = 0.05; n_steps = 4;   % on this coarse pair the fixed-Sigma march grows already for a static offset of 0.02
d_t = sprintf("%.17g*sin(2*pi*t)", A);
g = "(1 + t + t^2/2)";
us = ["(x - " + d_t + ")*(y - 0.5)*" + g; "((x - " + d_t + ")^2 - y)*" + g];
uu = "x*y*(1 + t)";
pp = "(x + y)*(1 + t)";
out = march(@(t) translating_circle_problems(k, curve_0, A, t, dt, uu, pp, us, ...
    rho_f, mu, rho_s, E, nu), fluid_mesh, solid_mesh, k, dt, n_steps, @(t) A*sin(2*pi*t));
fprintf("3. A = %g: worst L2 error per step %s, paths [%.3f, %.3f]\n", A, ...
    mat2str(out.errors, 3), out.paths);
n_fail = n_fail + check(max(out.errors) <= 1e-9, "the translating march is exact on the polynomial solution");
n_fail = n_fail + check(out.paths(1) < gap - A/2 && out.paths(2) > gap + A/2, ...
    "the transfer paths follow the moving interface");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_translating_interface:failed", "%d check(s) failed", n_fail);
end

function out = march(data, fluid_mesh, solid_mesh, k, dt, n_steps, d)
% Newmark solid + backward Euler, then Crank-Nicolson fluid, from
% FluidStructureCoupling.initial_state; d = [] keeps the history fixed,
% otherwise the history is translated by d(t_n) - d(t_n-1)
[fp, sp, src, ~, curve] = data(0);
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "euler";
state = FluidStructureCoupling.initial_state(fp, sp, fluid_mesh, solid_mesh, curve, src{:}, Order = k);
errors = zeros(1, n_steps);
paths = [Inf, 0];
for n = 1:n_steps
    t = n*dt;
    [fp, sp, src, pd, curve] = data(t);
    sp.Data.time_scheme = "newmark";
    fp.Data.time_scheme = "cranknicolson";
    if n == 1, fp.Data.time_scheme = "euler"; end
    if ~isempty(d)
        state = translate_state(sp, solid_mesh, state, [d(t) - d(t - dt); 0], k);
    end
    [fs, ss, state, gs] = FluidStructureCoupling.solve(fp, sp, fluid_mesh, solid_mesh, ...
        curve, src{:}, "newmark", quiet = true, state = state, patch_weight = "jacobian");
    e_F = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
    e_S = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
    errors(n) = max([e_F(:); e_S(:)]);
    for part = {gs.partition_1, gs.partition_2}
        l = vecnorm(part{1}.quadrature.interface - part{1}.quadrature.physical, 2, 1);
        paths = [min(paths(1), min(l(:))), max(paths(2), max(l(:)))];
    end
end
dofs = [{fs.Unknowns.field.dof}, {fs.Unknowns.hybrid.dof}, ...
        {ss.Unknowns.field.dof}, {ss.Unknowns.hybrid.dof}];
dofs = cellfun(@(c) double(cat(1, c{:})), cellfun(@(c) cellify(c), dofs, ...
    UniformOutput = false), UniformOutput = false);
out = struct("errors", errors, "paths", paths, "dofs", {dofs});
end

function c = cellify(x)
if iscell(x), c = x(:); else, c = {x}; end
end

function [fp, sp, src, pd, curve] = fixed_data(k, curve, t, dt, uu, pp, us, rho_f, mu, rho_s, E, nu)
% the fixed-Sigma data, load_problems on Sigma_0
[fp, sp, src, pd] = FluidStructureCoupling.load_problems(k, curve, t, dt, uu, pp, us, ...
    rho_f, mu, rho_s, E, nu, 0);
end

function s = project(solid, mesh, k, f)
% the volume and face projections of f and their lifted stress and rotation
s.displacement = L2proj2D(mesh, f, k);
s.displacement_skeleton = L2projSkeleton1d(mesh, f, k);
[s.stress, s.rotation] = solid.lift_stress(mesh, k, s.displacement, s.displacement_skeleton);
end
