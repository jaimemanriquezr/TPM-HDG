% test_fsi_hyperelastic_solid_path.m
% The solid side of the fsi kinematic row (InterfaceProblem.solve) for a
% HyperelasticProblem solid on a curved Sigma_0. The row extends the solid
% displacement from the face to Sigma_0 along the TPM path; for the
% hyperelastic solid the path gradient is its discrete displacement
% gradient L_h (HyperelasticProblem.discrete_gradient_operators), whereas
% the linear compliance c1 sigma - c2 tr(sigma) I + gamma misses
% -zeta^{-1} N(L_h): for a rigid rotation R, P = 0 and the linear
% reconstruction gives only the skew part of R - I, a defect (cos(theta) -
% 1) m along the path m.
% Exact rigid rotation u = (R - I)(X - c) of the disc of radius 0.4 about
% its centre c, theta = 0.5, Saint Venant-Kirchhoff, fluid at rest (exact
% zero), one shared circle, fitted_split meshes (chords: solid paths of
% length ~ h^2/(8R)), coupled Newton from the zero iterate:
%   0. the defect of the linear path reconstruction (path_gradient_integral
%      with sigma = 0, gamma = skew(R - I)) is (1 - cos(theta)) |m|;
%   1. steady (one backward-Euler step from rest, dt = 1), k = 1, 2: the
%      converged coupled solution is the exact one to round-off;
%   2. Newmark step (solid and fluid at density 0) from the converged
%      rotation theta_0 = 0.3 to theta = 0.5: the history of the
%      kinematic row (Data.state, transported by L_h) keeps the step exact,
%      and without Data.state the solve refuses the history;
%   3. smoke: lagrangian_carved_march on the circle (k = 1, A = 0.3, three
%      steps): every coupled Newton loop converges in at most 4 solves.
% With the linear compliance on the path the coupled errors of 1. are
% about 1e-3 at h = 1/8 and scale like the path length (h^2), for k = 1
% and 2 alike.
clc; clear;

n_fail = 0;
circle = CircularCurve([0.5; 0.5], 0.4);
c = [0.5; 0.5];
E = 10;  nu = 0.3;
theta = 0.5;

%% 0. defect of the linear path reconstruction
k = 2;
[fm, sm] = fitted_split(circle, 1/8);
sm = sm.add_stabilization(4/sqrt(2*max(sm.volume)));
[fp, sp, src] = FluidStructureCoupling.load_problems(k, circle, 0, 1, "0*x", "0*x", ...
    ["0*x"; "0*x"], 0, 1, 0, E, nu, 0);
[~, ~, ~, gs] = FluidStructureCoupling.solve(fp, sp, fm, sm, circle, src{:});
part = gs.partition_2;
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
c1 = 1/(2*mu);  c2 = lambda/(2*mu*(2*lambda + 2*mu));
G = get_affine_map(sm, part.elements, inverse = true);
gamma_dof = skew_dofs(k, sin(theta));
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
ratio = 0;
for e = 1:size(part.faces, 3)
    m = part.quadrature.interface(:, :, e) - part.quadrature.physical(:, :, e);
    if max(vecnorm(m)) < 1e-12, continue; end
    linear = path_gradient_integral(k, part, e, G(:, :, e), zeros(PEERSBasis(k).Dimension, 1), ...
        gamma_dof, c1, c2);
    ratio = max(ratio, max(vecnorm(linear.' - (R - eye(2))*m) ./ vecnorm(m)));
end
fprintf("linear path reconstruction: max |defect|/|m| = %.6f, 1 - cos(theta) = %.6f\n", ...
    ratio, 1 - cos(theta));
n_fail = n_fail + check(abs(ratio - (1 - cos(theta))) < 1e-10, ...
    "the linear reconstruction's defect is not (1 - cos(theta)) |m|");

%% 1. steady coupled rigid rotation
for k = [1, 2]
    [fp, sp, T1, T2, exact] = rotation_problems(k, circle, c, E, nu, theta, 0, 1);
    [fs, ss] = coupled_newton(fp, sp, fm, sm, circle, T1, T2, "steady", struct.empty);
    eS = ss.compute_l2_error(exact.solid);
    eF = fs.compute_l2_error(exact.fluid);
    fprintf("steady, k = %d: |u_S - u| %.2e, |gamma_S - gamma| %.2e, |P_h| %.2e, |u_F| %.2e\n", ...
        k, eS(2), eS(3), eS(1), eF(2));
    n_fail = n_fail + check(max([eS(:); eF(:)]) < 1e-10, ...
        sprintf("steady rigid rotation, k = %d: not exact", k));
end

%% 2. Newmark step with a history
k = 2;  dt = 0.1;  theta_0 = 0.3;
[fp, sp, T1, T2] = rotation_problems(k, circle, c, E, nu, theta_0, 0, 1);
[~, ss0] = coupled_newton(fp, sp, fm, sm, circle, T1, T2, "steady", struct.empty);
solid = struct("displacement", as_matrix(ss0.Unknowns.field(2).dof), ...
    "displacement_skeleton", as_matrix(ss0.Unknowns.hybrid(1).dof), ...
    "stress", as_matrix(ss0.Unknowns.field(1).dof), "rotation", as_matrix(ss0.Unknowns.field(3).dof));
state = struct("time", 0, "step", 0, "solid", solid, "fluid", struct.empty, "diagnostics", struct.empty);
[fp, sp, T1, T2, exact] = rotation_problems(k, circle, c, E, nu, theta, theta_0, dt);
[fs, ss] = coupled_newton(fp, sp, fm, sm, circle, T1, T2, "newmark", state);
eS = ss.compute_l2_error(exact.solid);
eF = fs.compute_l2_error(exact.fluid);
fprintf("Newmark step %.1f -> %.1f, k = %d: |u_S - u| %.2e, |gamma_S - gamma| %.2e, |P_h| %.2e, |u_F| %.2e\n", ...
    theta_0, theta, k, eS(2), eS(3), eS(1), eF(2));
n_fail = n_fail + check(max([eS(:); eF(:)]) < 1e-10, "Newmark rigid rotation step: not exact");
% the history needs the solid state: without Data.state the solve refuses
sp.Data.iterate = struct("displacement", solid.displacement, ...
    "displacement_skeleton", solid.displacement_skeleton);
n_fail = n_fail + check(throws(@() FluidStructureCoupling.solve(fp, sp, fm, sm, circle, T1, T2, ...
    "newmark", quiet = true, state = state), "InterfaceProblem:hyperelasticHistory"), ...
    "a hyperelastic history without Data.state was accepted");

%% 3. circle march smoke (lagrangian_carved_march passes the state)
background = StructuredSimplicialGrid(16, 16, "crossed");
r = lagrangian_carved_march(1, "circle", background, sm, 0.3, 0.1, 3, snap = true, ...
    solid_model = "hyperelastic");
fprintf("circle march, k = 1, A = 0.3: Newton iterations %s, |u_S - u| %.2e, |u_F - u| %.2e\n", ...
    mat2str([r.newton.iterations]), r.eS(2), r.eF(2));
n_fail = n_fail + check(all([r.newton.converged]) && max([r.newton.iterations]) <= 4 ...
    && r.eS(2) < 1e-3 && r.eF(2) < 1e-4, "circle hyperelastic march: Newton or errors off");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_hyperelastic_solid_path:failed", "%d check(s) failed", n_fail);
end

function [fp, sp, T1, T2, exact] = rotation_problems(k, circle, c, E, nu, theta, theta_0, dt)
% fluid at rest, hyperelastic solid rotated by theta about c (P = 0, no
% body force); T2 = u_F - (u_S - u_S0)/cdt with u_S0 the rotation by
% theta_0, cdt = dt (steady, theta_0 = 0) or dt/2 (trapezoidal Newmark)
[fp, ~, ~, pd] = FluidStructureCoupling.load_problems(k, circle, 0, dt, "0*x", "0*x", ...
    ["0*x"; "0*x"], 0, 1, 0, E, nu, 0);
sp = HyperelasticProblem(k, young_modulus = E, poisson_ratio = nu, time_step = dt);
rotation = @(a) [cos(a) - 1, -sin(a); sin(a), cos(a) - 1];
zero = @(x, y) 0*x;
D = rotation(theta);
cdt = dt;
if theta_0 ~= 0, cdt = dt/2; end
D_step = (D - rotation(theta_0))/cdt;
T1 = {zero; zero};
T2 = {@(x, y) -(D_step(1, 1)*(x - c(1)) + D_step(1, 2)*(y - c(2)));
      @(x, y) -(D_step(2, 1)*(x - c(1)) + D_step(2, 2)*(y - c(2)))};
u = {@(x, y) D(1, 1)*(x - c(1)) + D(1, 2)*(y - c(2)); @(x, y) D(2, 1)*(x - c(1)) + D(2, 2)*(y - c(2))};
% the rotation field: the (1,2) entry of skew(D) = [0 -s; s 0]
s = sin(theta);
gamma = @(x, y) -s + 0*x;
exact = struct("solid", {{repmat({zero}, 2, 2), u, gamma}}, ...
    "fluid", {{pd.sigma_F, pd.u_F, pd.gamma_F}});
end

function [fs, ss] = coupled_newton(fp, sp, fm, sm, circle, T1, T2, scheme, state)
% Newton on the coupled system from the zero iterate (as the march's
% coupled_newton, without damping): sp.Data.iterate <- the solved solid
k = sp.Order;
d1 = k + 1;  d2 = nchoosek(k + 2, 2);
sp.Data.iterate = struct("displacement", zeros(2*d2, size(sm.elements, 2)), ...
    "displacement_skeleton", zeros(2*d1, size(sm.faces, 2)));
if ~isempty(state)
    sp.Data.state = state;
    sp.Data.iterate = struct("displacement", state.solid.displacement, ...
        "displacement_skeleton", state.solid.displacement_skeleton);
end
for iteration = 1:15
    if isempty(state)
        [fs, ss] = FluidStructureCoupling.solve(fp, sp, fm, sm, circle, T1, T2, scheme, quiet = true);
    else
        [fs, ss] = FluidStructureCoupling.solve(fp, sp, fm, sm, circle, T1, T2, scheme, ...
            quiet = true, state = state);
    end
    u = as_matrix(ss.Unknowns.field(2).dof);
    change = max(abs(u - sp.Data.iterate.displacement), [], "all");
    sp.Data.iterate = struct("displacement", u, ...
        "displacement_skeleton", as_matrix(ss.Unknowns.hybrid(1).dof));
    if change < 1e-13, break; end
end
fprintf("  %s Newton: %d solves, last change %.1e\n", scheme, iteration, change);
end

function g = skew_dofs(k, s)
% Dubiner coefficients of the rotation field whose tensor is [0 -s; s 0]
A = asym_tensor_dubiner_basis(k, 1/3, 1/3);          % 1 x d2 x 2 x 2
g = zeros(nchoosek(k + 2, 2), 1);
g(1) = -s/A(1, 1, 1, 2);
end

function A = as_matrix(A)
if iscell(A), A = cat(1, A{:}); end
end

function failed = check(condition, message)
failed = ~condition;
if failed
    fprintf("FAILED: %s\n", message);
end
end
