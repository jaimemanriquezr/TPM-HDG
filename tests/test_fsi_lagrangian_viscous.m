% test_fsi_lagrangian_viscous.m
% The FSI march with a viscous hyperelastic Lagrangian solid
% (lagrangian_carved_march, solid_model = "hyperelastic", material =
% "neo_hookean", solid_viscosity = mu_v: the Kelvin-Voigt stress
% F mu_v (F_dot^T F + F^T F_dot) of HyperelasticProblem, assembled in the
% coupling by newton_blocks with the step's state): flat clamped Sigma_0,
% fluid carved from a "left" background grid and snapped, k = 1, N = 8,
% A = 0.3, T = 0.5, data from lagrangian_march_problems with the viscous
% stress in sigma_S and f_S (E = 10, so mu = 3.85; mu_v = 1).
%   1. every step converges and Newton is quadratic (consecutive relative
%      residuals with r_n < 0.2 and r_{n+1} > 1e-8: log r_{n+1}/log r_n
%      >= 1.7);
%   2. the errors against the exact solution at dt = 0.025 are at the
%      spatial level of this mesh (as without viscosity);
%   3. the term is significant: the solid stresses with mu_v = 1 and
%      mu_v = 0, each on its own data (same exact u_S, stresses differing
%      by F S_v), differ by more than 10x the mu_v = 1 stress error; with
%      2., the coupled scheme converges to the viscous solution.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
k = 1;  N = 8;  T = 0.5;  mu_v = 1;
background = StructuredSimplicialGrid(N, N, "left");
solid_mesh = apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
    "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
march = @(dt, viscosity) lagrangian_carved_march(k, "flat", background, solid_mesh, 0.3, dt, ...
    round(T/dt), "snap", true, solid_model = "hyperelastic", material = "neo_hookean", ...
    solid_viscosity = viscosity);

%% 1. Newton
r = march(0.1, mu_v);
orders = [];
for n = 1:numel(r.newton)
    h = r.newton(n).residual;
    fprintf("step %d (%s): %s\n", n, r.newton(n).stop_reason, sprintf("%.2e ", h));
    pairs = find(h(1:end-1) < 0.2 & h(2:end) > 1e-8);
    orders = [orders, log(h(pairs + 1))./log(h(pairs))]; %#ok<AGROW>
end
fprintf("observed orders log r_{n+1}/log r_n: %s\n", sprintf("%.2f ", orders));
n_fail = n_fail + check(all([r.newton.converged]), "a step did not converge");
n_fail = n_fail + check(~isempty(orders) && all(orders >= 1.7), "Newton not quadratic");

%% 2. errors at the spatial level
viscous = march(0.025, mu_v);
fprintf("mu_v = %g, dt = 0.025: |u_F - u| %.2e, |sigma_S - P| %.2e, |u_S - u| %.2e, geometry %.2e\n", ...
    mu_v, viscous.eF(2), viscous.eS(1), viscous.eS(2), max(viscous.history(:, 2)));
n_fail = n_fail + check(viscous.eF(2) < 1e-3 && viscous.eS(2) < 1e-3 && max(viscous.history(:, 2)) < 1e-3, ...
    "errors above the spatial level of N = 8");

%% 3. the viscous term is active
elastic = march(0.025, 0);
difference = max(abs(viscous.ss.Unknowns.field(1).dof - elastic.ss.Unknowns.field(1).dof), [], "all");
fprintf("max |sigma_S(mu_v = %g) - sigma_S(mu_v = 0)| = %.2e (error %.2e)\n", mu_v, difference, viscous.eS(1));
n_fail = n_fail + check(difference > 10*viscous.eS(1), "the viscous stress is not significant");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_viscous:failed", "%d check(s) failed", n_fail);
end
