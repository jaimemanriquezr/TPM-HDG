% test_fsi_eulerian_feedback_march.m
% eulerian_solid_march with interface_feedback = true (tests/helpers): both
% meshes carved and snapped against the interface predicted from the solid
% state (default predictor "none",carried interface traces), T1 on
% that curve (eulerian_interface_traction), the solved curve measured in
% out.feedback. Flat Sigma(t) = {y = 0.45 + d(t)}, d = A sin(2 pi t),
% A = 0.1, "left" grid.
%   1. T1 on a curve: on the exact Sigma(t) (the factory's curve and the
%      same set as a DeformedCurve over Sigma_0, flat and circle) the
%      handles equal load_problems' T_1 to 1e-12.
%   2. exactness: k = 2, N = 8, 6 steps of dt = 0.05, the polynomial data
%      of check 2 of test_fsi_eulerian_interface_march: worst L2 error per
%      step below 1e-8 and the solved curve's point error below 1e-8 at
%      every t_n, while the predicted curve lags the solved one.
%   3. smooth tangential data (eulerian_solid_convergence's solid_data =
%      "tangential"), k = 1, N = 8 and 16, T = 0.25, dt = 0.05 (8/N): the
%      march is bounded, eF and eS (sigma, u, gamma) within a factor 2 of
%      the prescribed march's, the lag (max over the steps) decreases with
%      N, the solved curve's point error within a factor 2 of the
%      measure-only one (out.interface of the prescribed march), and the
%      solved curve is monotone.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
A = 0.1;
d_t = sprintf("%.17g*sin(2*pi*t)", A);
uu = "sin(pi*x)*sin(pi*y)*cos(t)";  pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
us = ["0.1*sin(pi*x)*cos(y)*cos(t)"; "0.2*(y - 0.45 - " + d_t + ")*cos(x - y)*(1 + t^2)"];
physics = {1, 1, 1, 10, 0.3};

%% 1. T1 on the exact curve
t = 0.13;
c0 = CircularCurve([0.5; 0.5], 0.2);
cases = {@() translating_flat_problems(1, 0.45, A, t, 0.05, uu, pp, us, physics{:}), ...
    FlatCurveHorizontal(0.45), [0; 1]; ...
    @() translating_circle_problems(1, c0, A, t, 0.05, uu, pp, us, physics{:}), c0, [1; 0]};
names = ["flat", "circle"];
for c = 1:2
    [~, ~, src, pd, curve] = cases{c, 1}();
    shift = A*sin(2*pi*t)*cases{c, 3};
    deformed = DeformedCurve(cases{c, 2}, @(s) shift + 0*s, parametric = true, ...
        derivative = @(s) zeros(2, numel(s)));
    s = reshape(linspace(0.01, 0.99, 12), 3, 4);   % a 3 x 4 page of points
    gap = 0;
    for on = {curve, deformed}
        T1 = eulerian_interface_traction(on{1}, pd.sigma_F, pd.sigma_S);
        P = reshape(curve.get_2d_points(reshape(s, 1, [])), 2, []);
        x = reshape(P(1, :), size(s));  y = reshape(P(2, :), size(s));
        for i = 1:2
            gap = max(gap, max(abs(T1{i}(x, y) - src{1}{i}(x, y)), [], "all"));
        end
    end
    fprintf("1. %s: max |T1 - T_1| on the exact curve %.2e\n", names(c), gap);
    n_fail = n_fail + check(gap <= 1e-12, names(c) + " T1 differs from load_problems' T_1");
end

%% 2. exact polynomial solution with feedback
k = 2;  dt = 0.05;  n_steps = 6;
g = "(1 + t + t^2/2)";
us_poly = ["0.2*x*(1 - x)*" + g; "(y - 0.45 - " + d_t + ")*(0.3 + 0.2*x)*" + g];
background = StructuredSimplicialGrid(8, 8, "left");
[~, out] = evalc("eulerian_solid_march(k, background, A, dt, n_steps, ""x*y*(1 + t)"", " + ...
    """(x + y)*(1 + t)"", us_poly, verbose = true, interface_feedback = true)");
F = out.feedback;
fprintf("2. polynomial, k = 2: worst L2 error per step %s\n", mat2str(out.errors.', 3));
fprintf("   solved curve: max distance %.2e, max point error %.2e; lag max %.2e (normal %.2e)\n", ...
    max(F(:, 6)), max(F(:, 7)), max(F(2:end, 2)), max(F(2:end, 3)));
n_fail = n_fail + check(max(out.errors) <= 1e-8, "polynomial march with feedback not exact");
n_fail = n_fail + check(max(F(:, 6:7), [], "all") <= 1e-8, "polynomial solved curve not exact");
n_fail = n_fail + check(max(F(2:end, 2)) > 1e-5, "the predicted curve does not lag (test is void)");

%% 3. smooth tangential data
k = 1;  levels = [8 16];
lag = zeros(1, numel(levels));
for i = 1:numel(levels)
    N = levels(i);
    dt = 0.05*(levels(1)/N)^((k + 1)/2);
    n_steps = round(0.25/dt);
    background = StructuredSimplicialGrid(N, N, "left");
    prescribed = eulerian_solid_march(k, background, A, dt, n_steps, uu, pp, us, measure_interface = true);
    out = eulerian_solid_march(k, background, A, dt, n_steps, uu, pp, us, interface_feedback = true);
    F = out.feedback;
    ratio = [out.eF(1:3)./prescribed.eF(1:3), out.eS(1:3)./prescribed.eS(1:3)];
    point = max(F(2:end, 7))/max(prescribed.interface(2:end, 3));
    lag(i) = max(F(2:end, 2));
    fprintf("3. N = %2d: eF %s eS %s; ratios to prescribed %s\n", N, mat2str(out.eF(1:3), 3), ...
        mat2str(out.eS(1:3), 3), mat2str(ratio, 3));
    fprintf("   lag max %.3e (normal %.3e); solved point error max %.3e (ratio to measure-only %.2f), " + ...
        "min monotonicity %.2e, max|u_h| %.2f\n", lag(i), max(F(2:end, 3)), max(F(2:end, 7)), point, ...
        min(F(:, 9)), max(out.history(:, 2:3), [], "all"));
    n_fail = n_fail + check(all(isfinite(out.history(:, 2:3)), "all") && ...
        max(out.history(:, 2:3), [], "all") < 10, "feedback march unbounded");
    n_fail = n_fail + check(all(ratio <= 2), "feedback errors above twice the prescribed ones");
    n_fail = n_fail + check(point <= 2, "solved point error above twice the measure-only one");
    n_fail = n_fail + check(min(F(:, 9)) > 0, "solved curve not monotone");
end
fprintf("   lag reduction factor %.2f\n", lag(1)/lag(2));
n_fail = n_fail + check(lag(2) < lag(1), "lag does not decrease with N");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_eulerian_feedback_march:failed", "%d check(s) failed", n_fail);
end
