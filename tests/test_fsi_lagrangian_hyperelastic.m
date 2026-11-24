% test_fsi_lagrangian_hyperelastic.m
% The FSI march with a hyperelastic Lagrangian solid (lagrangian_carved_march,
% solid_model = "hyperelastic": the coupled Newton loop of each step, the
% solid's Newton step assembled in the coupling by
% HyperelasticProblem.newton_blocks): flat clamped Sigma_0, fluid carved
% from a "left" background grid and snapped, k = 1, N = 8, T = 0.5, data
% from lagrangian_march_problems (t^4 factor, sigma_S = P(I + grad u_S)).
%   1. quadratic Newton convergence in every step, A = 0.3 (|grad u_S| up
%      to about 0.2 at T), Saint Venant-Kirchhoff and neo-Hookean: every
%      step stops converged, and for consecutive relative residuals with
%      r_n < 0.2 and r_{n+1} > 1e-8 (above the round-off floor of the
%      coupled solves), log r_{n+1} / log r_n >= 1.7;
%   2. small amplitude: the hyperelastic march (SVK) against the linear
%      one on the same solid mesh (tau = max(30 mu, 2.4 lambda)/h on both,
%      the march's hyperelastic default), A = 0.02, 0.01,
%      0.005 with the fluid data scaled by A/0.02 (the whole solution
%      proportional to A): the solid dof difference scales as A^2 (log2
%      ratio in [1.8, 2.2]);
%   3. temporal self-convergence of order 2 of the solid (u, sigma) and of
%      the curve, A = 0.3, dt = 0.1, 0.05, 0.025 (bound 1.8 for the solid;
%      1.7 for the curve, whose velocity predictor error, of order dt^2,
%      is still pre-asymptotic at these steps);
%   4. the errors against the exact solution at dt = 0.025 are at the
%      spatial level of this mesh.
% Spatial convergence tables: experiments/moving-interface/
% lagrangian_march_convergence.m with solid_model = "hyperelastic".
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
k = 1;  N = 8;  T = 0.5;
background = StructuredSimplicialGrid(N, N, "left");
solid_mesh = apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
    "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
march = @(A, dt, varargin) lagrangian_carved_march(k, "flat", background, solid_mesh, A, dt, ...
    round(T/dt), "snap", true, varargin{:});

%% 1. quadratic Newton convergence per step
for material = ["svk", "neo_hookean"]
    r = march(0.3, 0.1, solid_model = "hyperelastic", material = material);
    orders = [];
    for n = 1:numel(r.newton)
        h = r.newton(n).residual;
        fprintf("%s step %d (%s): %s\n", material, n, r.newton(n).stop_reason, sprintf("%.2e ", h));
        pairs = find(h(1:end-1) < 0.2 & h(2:end) > 1e-8);
        orders = [orders, log(h(pairs + 1))./log(h(pairs))]; %#ok<AGROW>
    end
    fprintf("%s: observed orders log r_{n+1}/log r_n: %s\n", material, sprintf("%.2f ", orders));
    n_fail = n_fail + check(all([r.newton.converged]), material + ": a step did not converge");
    n_fail = n_fail + check(~isempty(orders) && all(orders >= 1.7), material + ": Newton not quadratic");
end

%% 2. small amplitude: hyperelastic - linear = O(A^2)
h = sqrt(2*max(solid_mesh.volume));
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(10, 0.3);   % lagrangian_march_problems
linear_mesh = solid_mesh.add_stabilization(max(30*mu, 2.4*lambda)/h);
amplitudes = [0.02, 0.01, 0.005];
difference = zeros(size(amplitudes));
for i = 1:numel(amplitudes)
    scale = amplitudes(i)/amplitudes(1);
    nonlinear = march(amplitudes(i), 0.1, solid_model = "hyperelastic", fluid_scale = scale);
    linear = lagrangian_carved_march(k, "flat", background, linear_mesh, amplitudes(i), 0.1, ...
        round(T/0.1), snap = true, fluid_scale = scale);
    difference(i) = max(abs(nonlinear.ss.Unknowns.field(2).dof - linear.ss.Unknowns.field(2).dof), [], "all");
end
ratio = log2(difference(1:end-1)./difference(2:end));
fprintf("max |u_S hyperelastic - u_S linear|: %s; log2 ratios %s\n", ...
    sprintf("%.3e ", difference), sprintf("%.2f ", ratio));
n_fail = n_fail + check(all(ratio > 1.8 & ratio < 2.2), "hyperelastic - linear not O(A^2)");

%% 3. temporal self-convergence
dts = [0.1, 0.05, 0.025];
runs = cell(1, 3);
for i = 1:3
    runs{i} = march(0.3, dts(i), solid_model = "hyperelastic");
end
s = linspace(0, 1, 401);
D = zeros(2, 3);
for i = 1:2
    a = runs{i};  b = runs{i + 1};
    D(i, :) = [max(abs(a.ss.Unknowns.field(2).dof - b.ss.Unknowns.field(2).dof), [], "all"), ...
               max(abs(a.ss.Unknowns.field(1).dof - b.ss.Unknowns.field(1).dof), [], "all"), ...
               max(vecnorm(a.curve.get_2d_points(s) - b.curve.get_2d_points(s)))];
end
order = log2(D(1, :) ./ D(2, :));
fprintf("self-convergence order (u_S, sigma_S, Sigma_h): %.2f %.2f %.2f\n", order);
n_fail = n_fail + check(all(order(1:2) > 1.8) && order(3) > 1.7, "temporal self-convergence below order 2");

%% 4. errors at the spatial level
r = runs{3};
fprintf("dt = %g: |u_F - u| %.2e, |sigma_S - P| %.2e, |u_S - u| %.2e, geometry %.2e\n", ...
    dts(3), r.eF(2), r.eS(1), r.eS(2), max(r.history(:, 2)));
n_fail = n_fail + check(r.eF(2) < 1e-3 && r.eS(2) < 1e-3 && max(r.history(:, 2)) < 1e-3, ...
    "errors above the spatial level of N = 8");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_hyperelastic:failed", "%d check(s) failed", n_fail);
end
