% test_fsi_eulerian_interface_march.m
% eulerian_solid_march with measure_interface = true (tests/helpers): after
% every step the interface is rebuilt from the solved solid displacement
% (eulerian_interface_displacement from the P_k+1 postprocessing u* and
% from u_h) and measured against the exact set Sigma_0 + d(t) and the
% exact material points (interface_curve_metrics); nothing is fed back.
% Flat Sigma(t) = {y = 0.45 + d(t)}, d = A sin(2 pi t), A = 0.1, "left"
% grid N = 8, snapped. The solid data map Sigma(t) into Sigma_0 (zero
% normal component on Sigma(t), zero tangential one at the walls x = 0, 1),
% so the exact interface is the prescribed one while the material slides
% along it.
%   1. measuring does not touch the solve: eF, eS equal (isequal) those of
%      the march without it (k = 1, 3 steps of dt = 0.05, the smooth data of check 3).
%   2. exactness: k = 2, 6 steps of dt = 0.05 (the setting of check 3 of
%      test_fsi_eulerian_solid_march, its fluid data), u_S = [0.2 x (1 - x);
%      (y - 0.45 - d)(0.3 + 0.2 x)] g(t), g = 1 + t + t^2/2: worst L2
%      error per step below 1e-9, and the curve rebuilt from u* has
%      distance and point error below 1e-9 at every t_n (t = 0 included),
%      while the material slides tangentially by more than 0.03.
%   3. smooth data, u_S = [0.1 sin(pi x) cos(y) cos(t); 0.2 (y - 0.45 -
%      d) cos(x - y)(1 + t^2)], k = 1, N = 8 and 16, T = 0.25, dt = 0.05
%      (8/N) (the driver's dt0 (N0/N)^((k+1)/2)): the maxima over the
%      steps of the distance and the point error of the u* curve decrease
%      by factors above 4 and 3 (measured 5.83 and 3.52; with
%      skeleton_history = "library" 7.56 and 5.23: the u* trace of the
%      default displacement skeleton history slows the point error), the
%      continuous curve is monotone and its measure ratio within 1e-3 of 1
%      at every step.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
A = 0.1;
background = StructuredSimplicialGrid(8, 8, "left");
d_t = sprintf("%.17g*sin(2*pi*t)", A);
uu = "sin(pi*x)*sin(pi*y)*cos(t)";  pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
us = ["0.1*sin(pi*x)*cos(y)*cos(t)"; "0.2*(y - 0.45 - " + d_t + ")*cos(x - y)*(1 + t^2)"];

%% 1. measuring does not change the solve
out = eulerian_solid_march(1, background, A, 0.05, 3, uu, pp, us);
measured = eulerian_solid_march(1, background, A, 0.05, 3, uu, pp, us, measure_interface = true);
same = isequal(out.eF, measured.eF) && isequal(out.eS, measured.eS);
fprintf("1. eF, eS with and without measure_interface identical: %d (%d rows measured)\n", ...
    same, size(measured.interface, 1));
n_fail = n_fail + check(same && size(measured.interface, 1) == 4, "measure_interface changes the solve");

%% 2. exact polynomial solution, exact interface
k = 2;  dt = 0.05;  n_steps = 6;
g = "(1 + t + t^2/2)";
us_poly = ["0.2*x*(1 - x)*" + g; "(y - 0.45 - " + d_t + ")*(0.3 + 0.2*x)*" + g];
[~, out] = evalc("eulerian_solid_march(k, background, A, dt, n_steps, ""x*y*(1 + t)"", " + ...
    """(x + y)*(1 + t)"", us_poly, verbose = true, measure_interface = true)");   % per-step errors, unprinted
worst = out.errors.';
iface = out.interface;
% the tangential slide at the final time: the exact material point x(s) of
% X(s) = (s, 0.45) minus X(s) - d e_y (u_S1 does not depend on y)
s = linspace(0, 1, 101);
T = n_steps*dt;
slide = max(abs(0.2*s.*(1 - s)*(1 + T + T^2/2)));
fprintf("2. polynomial, k = 2: worst L2 error per step %s\n", mat2str(worst, 3));
fprintf("   interface from u*: max distance %.2e, max point error %.2e over t = 0..%g; slide %.3f\n", ...
    max(iface(:, 2)), max(iface(:, 3)), T, slide);
n_fail = n_fail + check(max(worst) <= 1e-9, "polynomial march exact");
n_fail = n_fail + check(max(iface(:, 2:3), [], "all") <= 1e-9, "polynomial interface exact");
n_fail = n_fail + check(slide > 0.03, "polynomial interface slides by less than 0.03");

%% 3. convergence of the measured interface, smooth data
k = 1;  levels = [8 16];
E = zeros(numel(levels), 4);
for i = 1:numel(levels)
    N = levels(i);
    dt = 0.05*(levels(1)/N)^((k + 1)/2);
    n_steps = round(0.25/dt);
    out = eulerian_solid_march(k, StructuredSimplicialGrid(N, N, "left"), A, dt, n_steps, uu, pp, us, ...
        measure_interface = true);
    steps = out.interface(2:end, :);
    E(i, :) = [max(steps(:, 2)), max(steps(:, 3)), max(abs(out.interface(:, 4) - 1)), ...
        min(out.interface(:, 5))];
    fprintf("3. N = %2d dt %.4f: max distance %.3e, max point error %.3e, max |ratio - 1| %.2e, " + ...
        "min monotonicity %.2e (raw %.2e), max jump %.2e; u_h curve %.3e / %.3e\n", N, dt, E(i, :), ...
        min(out.interface(:, 6)), max(steps(:, 7)), max(steps(:, 8)), max(steps(:, 9)));
end
factor = E(1, 1:2)./E(2, 1:2);
fprintf("   reduction factors (distance, point error) %.2f %.2f\n", factor);
n_fail = n_fail + check(all(factor > [4, 3]), "interface distance / point error reduced by less than 4 / 3");
n_fail = n_fail + check(all(E(:, 4) > 0), "continuous interface not monotone");
n_fail = n_fail + check(all(E(:, 3) < 1e-3), "measure ratio off by 1e-3");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_eulerian_interface_march:failed", "%d check(s) failed", n_fail);
end
