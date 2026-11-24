% test_fsi_lagrangian_exact_geometry.m
% The Lagrangian-interface FSI march (lagrangian_carved_march) with the
% coupling data on the exact curve X + u_S (lagrangian_march_problems,
% exact_geometry = true), so that the errors include the march's geometry
% error (predictor, lag, reconstruction of Sigma_h). Flat Sigma_0, k = 1,
% N = 8, T = 0.5, dt = 0.1, 0.05, 0.025, snapped fluid.
%   1. exact_geometry = false (and predictor = "velocity") given explicitly
%      reproduces the default march bit for bit;
%   2. with exact_geometry = true the errors against the exact solution
%      stay within 3 times those of the data on the discrete curve;
%   3. temporal self-convergence of the solid (u, sigma) and of the curve
%      with exact_geometry = true, order > 1.6 for the Newmark predictor
%      (printed, not asserted, for the Taylor and velocity predictors), and the
%      difference of the solid between the two data (the error the
%      discrete-curve data hide);
%   4. geometry_iterations = 1 against 0 with exact_geometry = true
%      (printed).
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
k = 1;  N = 8;  T = 0.5;  A = 0.03;
background = StructuredSimplicialGrid(N, N, "left");
solid_mesh = apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
    "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
data = @(dt, on_exact) @(t) lagrangian_march_problems(k, "flat", t, dt, A, exact_geometry = on_exact);
march = @(dt, varargin) lagrangian_carved_march(k, "flat", background, solid_mesh, A, dt, ...
    round(T/dt), "snap", true, varargin{:});
dts = [0.1, 0.05, 0.025];
predictors = ["newmark", "taylor", "velocity"];
discrete = cell(1, 3);
exact = cell(3, 3);
for i = 1:3
    discrete{i} = march(dts(i));
    for p = 1:3
        exact{p, i} = march(dts(i), "problems", data(dts(i), true), "predictor", predictors(p));
    end
end

%% 1. the default data and predictor, given explicitly
explicit = march(dts(1), "problems", data(dts(1), false), "predictor", "velocity");
a = discrete{1};
same = isequal(a.eF, explicit.eF) && isequal(a.eS, explicit.eS) && isequal(a.history, explicit.history) ...
    && isequal(a.paths, explicit.paths) ...
    && isequal(a.ss.Unknowns.field(2).dof, explicit.ss.Unknowns.field(2).dof) ...
    && isequal(a.fs.Unknowns.field(2).dof, explicit.fs.Unknowns.field(2).dof);
fprintf("exact_geometry = false, predictor = ""velocity"" bit-identical to the default: %d\n", same);
n_fail = n_fail + check(same, "explicit default options change the march");

%% 2.-3. errors and self-convergence
s = linspace(0, 1, 401);
difference = @(a, b) [max(abs(a.ss.Unknowns.field(2).dof - b.ss.Unknowns.field(2).dof), [], "all"), ...
    max(abs(a.ss.Unknowns.field(1).dof - b.ss.Unknowns.field(1).dof), [], "all"), ...
    max(vecnorm(a.curve.get_2d_points(s) - b.curve.get_2d_points(s)))];
fprintf("\n%-9s %-8s %-6s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n", "data", "pred.", "dt", ...
    "|uF-u|", "|sF-s|", "|uS-u|", "|sS-s|", "geometry", "predictor", "|uS-uS_d|");
labels = ["discrete", "exact", "exact", "exact"];
series = {discrete, exact(1, :), exact(2, :), exact(3, :)};
names = ["velocity", "newmark", "taylor", "velocity"];
for j = 1:4
    for i = 1:3
        r = series{j}{i};
        hidden = difference(r, discrete{i});
        fprintf("%-9s %-8s %-6.3g %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e\n", ...
            labels(j), names(j), dts(i), r.eF(2), r.eF(1), r.eS(2), r.eS(1), ...
            max(r.history(:, 2)), max(r.history(:, 3)), hidden(1));
    end
end
fprintf("(geometry, predictor: max over steps of max |Sigma_h - (X + u_S)| on Sigma_0 for the solved\n" + ...
    " and the predicted curve; |uS-uS_d|: max solid u dof difference to the discrete-curve data)\n");

e = exact{1, 3};  d = discrete{3};
ratio = [e.eF(2)/d.eF(2), e.eS(2)/d.eS(2), max(e.history(:, 2))/max(d.history(:, 2))];
fprintf("\nexact/discrete data at dt = %g (uF, uS, geometry): %.3f %.3f %.3f\n", dts(3), ratio);
n_fail = n_fail + check(all(ratio < 3), "exact-geometry errors more than 3x the discrete-curve ones");

fprintf("\nself-convergence order, dt = 0.1/0.05/0.025 (u_S, sigma_S, Sigma_h):\n");
order = zeros(3, 3);
for p = 1:3
    D = [difference(exact{p, 1}, exact{p, 2}); difference(exact{p, 2}, exact{p, 3})];
    order(p, :) = log2(D(1, :) ./ D(2, :));
    fprintf("  %-8s %.2f %.2f %.2f\n", predictors(p), order(p, :));
end
H = [difference(exact{1, 1}, discrete{1}); difference(exact{1, 2}, discrete{2}); ...
     difference(exact{1, 3}, discrete{3})];
fprintf("hidden difference (newmark, exact - discrete data), u_S: %s; ratios %s\n", ...
    sprintf("%.2e ", H(:, 1)), sprintf("%.2f ", H(1:2, 1) ./ H(2:3, 1)));
n_fail = n_fail + check(all(order(1, :) > 1.6), "exact-geometry self-convergence below 1.6 (newmark)");

%% 4. one geometry iteration with exact_geometry
it = march(dts(1), "problems", data(dts(1), true), "geometry_iterations", 1);
r = exact{1, 1};
change = [it.eF(1:2), it.eS(1:2)] ./ [r.eF(1:2), r.eS(1:2)] - 1;
fprintf("\ngeometry_iterations 1 vs 0, exact data, dt = %g: lag %.2e -> %.2e; " + ...
    "relative error change (sF uF sS uS) %s; |uS(1) - uS(0)| %.2e\n", dts(1), ...
    max(r.history(:, 8)), max(it.history(:, 8)), sprintf("%+.1e ", change), difference(it, r)*[1; 0; 0]);

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_exact_geometry:failed", "%d check(s) failed", n_fail);
end
