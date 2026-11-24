% test_fsi_eulerian_carved_march.m
% eulerian_carved_march (tests/helpers): the fluid re-carved from a fixed
% background grid at every step, its history remapped by background
% element, the solid translating rigidly, k = 1.
%   1. flat, Sigma(t) = {y = 0.5 + 0.1 sin(2 pi t)}, "right" grid N = 8,
%      structured solid, t^4 data (translating_flat_problems), dt = 0.05,
%      T = 0.75 (the interface rises by most of an element, returns and
%      descends as far): elements switch on, by whole rows of the grid
%      (2N elements each) as the interface crosses a grid line, the
%      solution stays bounded and the errors are at the level of this
%      coarse mesh;
%   2. circle with snap = true (snap_outside), crossed grid N = 16, A =
%      0.02, dt = 0.05, 3 steps, the fitted_split disc (mesh2d): the march
%      completes with bounded errors, and the switched-on / off counts of
%      every step equal those of the unsnapped carves (the background
%      elements that enter / leave the fluid), below 10 % of the fluid
%      elements.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
k = 1;

%% 1. flat interface
N = 8;  A = 0.1;  dt = 0.05;  n_steps = 15;
uf = "sin(pi*x)*sin(pi*y)*t^4";  pf = "cos(pi*x)*(y - 0.5)*t^4";
sf = ["0.1*sin(2*x + y)*t^4"; "0.1*cos(x - y)*t^4"];
background = StructuredSimplicialGrid(N, N, "right");
solid_mesh = apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
    "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
timer = tic;
out = eulerian_carved_march(k, background, solid_mesh, A, dt, n_steps, uf, pf, sf, geometry = "flat");
fprintf("1. flat, N = %d, A = %g, T = %g (%.1f s): |uF - u| %.3e, |uS - u| %.3e, |sF - s| %.3e, |sS - s| %.3e\n", ...
    N, A, dt*n_steps, toc(timer), out.eF(2), out.eS(2), out.eF(1), out.eS(1));
fprintf("   max |u_h| (fluid, solid) %.3g %.3g; switched on %s, off %s; paths [%.2e, %.2e]\n", ...
    max(out.history(:, 1)), max(out.history(:, 2)), mat2str(out.history(:, 3).'), ...
    mat2str(out.history(:, 4).'), out.paths);
n_fail = n_fail + check(all(isfinite([out.eF, out.eS])) && max(out.history(:, 1:2), [], "all") < 10, ...
    "flat: the march is unbounded");
n_fail = n_fail + check(out.eF(2) < 3e-2 && out.eS(2) < 1e-3, "flat: errors above the level of N = 8");
n_fail = n_fail + check(sum(out.history(:, 3)) > 0, "flat: no element switched on");
n_fail = n_fail + check(all(mod(out.history(:, 3:4), 2*N) == 0, "all"), ...
    "flat: elements switched on / off other than by whole rows of the grid");

%% 2. circle, snapped
assert(exist("refine2", "file") == 2, "test_fsi_eulerian_carved_march:mesh2d", ...
    "mesh2d (refine2) is not on the path: add src/mesh-generation/utilities/mesh2d " + ...
    "(populated in the main checkout) before running this test.");
N = 16;  A = 0.02;  dt = 0.05;  n_steps = 3;
uu = "sin(pi*x)*sin(pi*y)*cos(t)";  pp = "cos(pi*x)*(y - 0.5)*(1 + t)";
us = ["0.1*sin(2*x + y)*cos(t)"; "0.1*cos(x - y)*(1 + t^2)"];
background = StructuredSimplicialGrid(N, N, "crossed");
[~, solid_mesh] = fitted_split(CircularCurve([0.5; 0.5], 0.4), 0.1);
timer = tic;
out = eulerian_carved_march(k, background, solid_mesh, A, dt, n_steps, uu, pp, us, snap = true);
n_fluid = size(out.fs.Mesh.elements, 2);
fprintf("2. circle snapped, N = %d, A = %g, %d steps of %g (%.1f s): |uF - u| %.3e, |uS - u| %.3e, |sF - s| %.3e, |sS - s| %.3e\n", ...
    N, A, n_steps, dt, toc(timer), out.eF(2), out.eS(2), out.eF(1), out.eS(1));
% the background elements entering / leaving the fluid: the unsnapped carves
centroid = @(m) squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
carve = @(t) background.split(CircularCurve([0.5 + A*sin(2*pi*t); 0.5], 0.4));
expected = zeros(n_steps, 2);
old = carve(0);
for n = 1:n_steps
    new = carve(n*dt);
    present = ismembertol(centroid(new), centroid(old), 1e-10, "ByRows", true);
    expected(n, :) = [nnz(~present), size(old.elements, 2) - nnz(present)];
    old = new;
end
fprintf("   switched on/off per step: march %s, background %s; %d fluid elements; max |u_h| %.3g %.3g\n", ...
    mat2str(out.history(:, 3:4)), mat2str(expected), n_fluid, max(out.history(:, 1)), max(out.history(:, 2)));
n_fail = n_fail + check(all(isfinite([out.eF, out.eS])) && max(out.history(:, 1:2), [], "all") < 10, ...
    "circle: the march is unbounded");
n_fail = n_fail + check(out.eF(2) < 1.8e-2 && out.eS(2) < 1.7e-3, "circle: errors above the level of N = 16");
n_fail = n_fail + check(isequal(out.history(:, 3:4), expected), ...
    "circle: switched counts differ from the background elements entering / leaving the fluid");
n_fail = n_fail + check(max(out.history(:, 3)) < 0.1*n_fluid, ...
    sprintf("circle: more than 10%% of the %d fluid elements switched on in a step", n_fluid));

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_eulerian_carved_march:failed", "%d check(s) failed", n_fail);
end
