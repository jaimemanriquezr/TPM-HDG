% Transient FSI with a prescribed translating interface on moving meshes
% (ale_translating_march): the solid mesh translates with Sigma(t) =
% Sigma_0 + A sin(2 pi t) e_x, the fluid mesh follows it by an x-ramp
% (ALE, beta = -w_h in the advective form), so the transfer paths keep
% their t = 0 lengths.
%
%   1. A = 0 reproduces the fixed-mesh march (load_problems on Sigma_0)
%      dof for dof (1e-12), at the solver's default k + 1 patch points.
%   2. A > 0 is exact on a solution the integrators and the ALE form
%      reproduce: u_F = (-2y(1 + t), 0) (x-independent, so its ALE
%      derivative is d_t u_F; linear in t), u_S = p(x - d(t)) g(t) with p of
%      degree k and g quadratic (Newmark at fixed material points). The
%      error stays below 1e-9 and the path lengths do not change.
%
%See also: ale_translating_march, test_fsi_translating_interface
clc; clear;

k = 2; h = 0.1; gap = 0.05;
rho_f = 1; mu = 1; rho_s = 1; E = 10; nu = 0.3;
curve_0 = CircularCurve([0.5; 0.5], 0.4);
[fluid_mesh, solid_mesh] = fitted_split(curve_0, h, gap);
n_fail = 0;

%% 1. A = 0 against the fixed-mesh march
uu = "sin(pi*x)*sin(pi*y)*cos(t)";
pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
us = ["x*sin(y)*cos(t)"; "y*cos(x)*(1 + t^2)"];
dt = 0.05; n_steps = 3;
ref = fixed_march(k, curve_0, fluid_mesh, solid_mesh, dt, n_steps, uu, pp, us, ...
    rho_f, mu, rho_s, E, nu);
ale = ale_translating_march(k, fluid_mesh, solid_mesh, 0, dt, n_steps, uu, pp, us, ...
    rho_f = rho_f, mu = mu, rho_s = rho_s, E = E, nu = nu, n_quadpoints = k + 1);
d_fix = max(cellfun(@(a, b) max(abs(a - b), [], "all"), dofs_of(ref.fs, ref.ss), ...
    dofs_of(ale.fs, ale.ss)));
fprintf("1. A = 0 vs fixed meshes: max dof difference %.2e\n", d_fix);
n_fail = n_fail + check(d_fix <= 1e-12, "A = 0 reproduces the fixed-mesh march");

%% 2. A > 0, exact solution
A = 0.02; dt = 0.05; n_steps = 4;
d_t = sprintf("%.17g*sin(2*pi*t)", A);
g = "(1 + t + t^2/2)";
us = ["(x - " + d_t + ")*(y - 0.5)*" + g; "((x - " + d_t + ")^2 - y)*" + g];
out = ale_translating_march(k, fluid_mesh, solid_mesh, A, dt, n_steps, ...
    "y^2*(1 + t)", "(x + y)*(1 + t)", us, rho_f = rho_f, mu = mu, rho_s = rho_s, E = E, nu = nu);
rest = ale_translating_march(k, fluid_mesh, solid_mesh, 0, dt, 1, ...
    "y^2*(1 + t)", "(x + y)*(1 + t)", us, rho_f = rho_f, mu = mu, rho_s = rho_s, E = E, nu = nu);
worst = max([out.eF, out.eS]);
fprintf("2. A = %g: worst L2 error %.2e, paths [%.4f, %.4f] (A = 0: [%.4f, %.4f])\n", ...
    A, worst, out.paths, rest.paths);
n_fail = n_fail + check(worst <= 1e-9, "the ALE march is exact on the polynomial solution");
n_fail = n_fail + check(max(abs(out.paths - rest.paths)) <= 1e-10, ...
    "the transfer paths keep their t = 0 lengths");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_ale_translating_interface:failed", "%d check(s) failed", n_fail);
end

function out = fixed_march(k, curve, fluid_mesh, solid_mesh, dt, n_steps, uu, pp, us, rho_f, mu, rho_s, E, nu)
% the fixed-mesh march: Newmark solid, backward Euler then Crank-Nicolson fluid
data = @(t) FluidStructureCoupling.load_problems(k, curve, t, dt, uu, pp, us, ...
    rho_f, mu, rho_s, E, nu, 0);
[fp, sp, src] = data(0);
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "euler";
state = FluidStructureCoupling.initial_state(fp, sp, fluid_mesh, solid_mesh, curve, src{:}, Order = k);
for n = 1:n_steps
    [fp, sp, src] = data(n*dt);
    sp.Data.time_scheme = "newmark";
    fp.Data.time_scheme = "cranknicolson";
    if n == 1, fp.Data.time_scheme = "euler"; end
    [fs, ss, state] = FluidStructureCoupling.solve(fp, sp, fluid_mesh, solid_mesh, curve, ...
        src{:}, "newmark", quiet = true, state = state, patch_weight = "jacobian");
end
out = struct("fs", fs, "ss", ss);
end

function dofs = dofs_of(fs, ss)
% every field and trace dof array of both solutions
fields = [{fs.Unknowns.field.dof}, {fs.Unknowns.hybrid.dof}, ...
          {ss.Unknowns.field.dof}, {ss.Unknowns.hybrid.dof}];
dofs = cellfun(@(c) double(cat(1, c{:})), cellfun(@cellify, fields, UniformOutput = false), ...
    UniformOutput = false);
end

function c = cellify(x)
if iscell(x), c = x(:); else, c = {x}; end
end
