% test_fsi_lagrangian_march.m
% The FSI march with the interface following the discrete Lagrangian solid
% (lagrangian_carved_march): flat Sigma_0, fluid carved from a "left"
% background grid and snapped, solid Newmark on its fixed reference mesh,
% data from lagrangian_march_problems (t^4 factor), k = 1, N = 8, T = 0.5.
%   1. temporal self-convergence of order 2 of the solid (u, sigma) and of
%      the curve, dt = 0.1, 0.05, 0.025 (the solid mesh is fixed, so its
%      dofs compare directly);
%   2. the march stays bounded and the errors against the exact solution
%      are at the spatial level of this mesh;
%   3. geometry_iterations = 1 cuts the lag of Sigma_h behind the solved
%      displacement by orders of magnitude and leaves the errors unchanged
%      to 1%.
% Convergence tables: experiments/moving-interface/lagrangian_march_convergence.m.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
k = 1;  N = 8;  T = 0.5;  A = 0.03;
background = StructuredSimplicialGrid(N, N, "left");
solid_mesh = apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
    "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
dts = [0.1, 0.05, 0.025];
runs = cell(1, 3);
for i = 1:3
    runs{i} = lagrangian_carved_march(k, "flat", background, solid_mesh, A, dts(i), ...
        round(T/dts(i)), snap = true);
end

%% 1. temporal self-convergence
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
n_fail = n_fail + check(all(order > 1.8), "temporal self-convergence below order 2");

%% 2. bounded, errors at the spatial level
r = runs{3};
fprintf("dt = %g: |u_F - u| %.2e, |u_S - u| %.2e, geometry %.2e, max |u_h| %.3g\n", ...
    dts(3), r.eF(2), r.eS(2), max(r.history(:, 2)), max(r.history(:, 4:5), [], "all"));
n_fail = n_fail + check(r.eF(2) < 1e-3 && r.eS(2) < 1e-3 && max(r.history(:, 2)) < 2e-3, ...
    "errors above the spatial level of N = 8");

%% 3. one geometry iteration
it = lagrangian_carved_march(k, "flat", background, solid_mesh, A, dts(1), round(T/dts(1)), ...
    snap = true, geometry_iterations = 1);
lag = [max(runs{1}.history(:, 8)), max(it.history(:, 8))];
change = max(abs([it.eF(1:2), it.eS(1:2)] ./ [runs{1}.eF(1:2), runs{1}.eS(1:2)] - 1));
fprintf("lag of Sigma_h: explicit %.2e, one iteration %.2e; relative error change %.1e\n", lag, change);
n_fail = n_fail + check(lag(2) < 1e-2*lag(1) && change < 1e-2, "geometry iteration");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_march:failed", "%d check(s) failed", n_fail);
end
